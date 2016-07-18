// this file is distributed under 
// MIT license
#ifndef ____table_data_H_____
#	define ____table_data_H_____
#include "sigma.h"
#include "points.h"
#include "chains.h"
#include "error.h"
namespace MathTemplates{
	template<class numX, class numY=numX>
	class SortedPoints:public SortedChain<point_editable_y<numX,numY>>{
	public:
		typedef std::function<numY(const numX&)> Func;
		SortedPoints(){}
		SortedPoints(const SortedChain<point<numX,numY>>&chain){
			for(const auto&x:chain)
				SortedChain<point_editable_y<numX,numY>>::append_item_from_sorted(x);
		}
		SortedPoints(const SortedChain<point<numX,numY>>&&chain):SortedPoints(chain){}
		SortedPoints(const Func f,const SortedChain<numX>&chain){
			for(numX x:chain)
				SortedChain<point_editable_y<numX,numY>>::append_item_from_sorted(point<numX,numY>(x,f(x)));
		}
		SortedPoints(const Func f,const SortedChain<numX>&&chain):SortedPoints(f,chain){}
		SortedPoints&operator<<(const point<numX,numY>&p){
			SortedChain<point_editable_y<numX,numY>>::operator<<(point_editable_y<numX,numY>(p));
			return *this;
		}
		SortedPoints&operator<<(const point<numX,numY>&&p){return operator<<(p);}
		SortedPoints(const std::initializer_list<point<numX,numY>>&points){
			for(const auto&p:points)operator<<(p);
		}
		SortedPoints(const Func f,const std::initializer_list<numX>&&chain){
			for(numX x:chain)operator<<(point<numX,numY>(x,f(x)));
		}
		
		virtual ~SortedPoints(){}
		
		SortedPoints(const SortedChain<point_editable_y<numX,numY>>&points):SortedChain<point_editable_y<numX,numY>>(points){}
		SortedPoints(const SortedChain<point_editable_y<numX,numY>>&&points):SortedPoints(points){}
		SortedPoints& operator=(const SortedChain<point_editable_y<numX,numY>>&points){
			SortedChain<point_editable_y<numX,numY>>::operator=(points);
			return *this;
		}
		
		point_editable_y<numX,numY>&Bin(const size_t i){return SortedChain<point_editable_y<numX,numY>>::accessBin(i);}

		const  std::vector<point<numY,numX>> Transponate()const{
			std::vector<point<numY,numX>> res;
			for(const auto&p:*this)
				res.push_back(point<numX,numY>(p.Y(),p.X()));
			return res;
		}
		const  SortedChain<point<numY,numX>> TransponateAndSort()const{
			SortedChain<point<numY,numX>> res;
			for(const auto&p:*this)
				res<<point<numX,numY>(p.Y(),p.X());
			return res;
		}
		const SortedPoints XRange(const numX from,const numX to)const{
			SortedPoints res;
			for(const auto&P:*this)
				if((P.X()>=from)&&(P.X()<=to))
					res<<P;
				return res;
		}
		const SortedPoints YRange(const numY from,const numY to)const{
			SortedPoints res;
			for(const auto&P:*this)
				if((P.Y()>=from)&&(P.Y()<=to))
					res<<P;
				return res;
		}
		const SortedPoints XExclude(const numX from,const numX to)const{
			SortedPoints res;
			for(const auto&P:*this)
				if((P.X()<from)||(P.X()>to))
					res<<P;
				return res;
		}
		const SortedPoints YExclude(const numY from,const numY to)const{
			SortedPoints res;
			for(const auto&P:*this)
				if((P.Y()<from)||(P.Y()>to))
					res<<P;
				return res;
		}
		SortedPoints&FillWithValues(const numY&v){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()=v;
			return *this;
		}
		SortedPoints&FillWithValues(const numY&&v){
			return FillWithValues(v);
		}
		SortedPoints&Transform(const std::function<numY(const numX&,const numY&)>&F){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()=F(Bin(i).X(),Bin(i).Y());
			return *this;
		}
		SortedPoints&operator+=(const SortedPoints& second){
			for(size_t i=0,n=this->size();i<n;i++){
				if(Bin(i).X()==second[i].X())
					Bin(i).varY()+=second[i].Y();
				else
					throw Exception<SortedPoints>("Cannot perform arithmetic actions two chains. X coordinates differ");
			}
			return *this;
		}
		SortedPoints&operator+=(const SortedPoints&&second){return operator+=(second);}
		SortedPoints&operator+=(const Func f){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()+=f(Bin(i).X());
			return *this;
		}
		SortedPoints&operator+=(const numY&c){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()+=c;
			return *this;
		}
		SortedPoints&operator+=(const numY&&c){
			return operator+=(c);
		}
		
		SortedPoints&operator-=(const SortedPoints& second){
			for(size_t i=0,n=this->size();i<n;i++){
				if(Bin(i).X()==second[i].X())
					Bin(i).varY()-=second[i].Y();
				else
					throw Exception<SortedPoints>("Cannot perform arithmetic actions two chains. X coordinates differ");
			}
			return *this;
		}
		SortedPoints&operator-=(const SortedPoints&&second){return operator-=(second);}
		SortedPoints&operator-=(const Func f){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()-=f(Bin(i).X());
			return *this;
		}
		SortedPoints&operator-=(const numY&c){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()-=c;
			return *this;
		}
		SortedPoints&operator-=(const numY&&c){
			return operator-=(c);
		}
		
		SortedPoints&operator*=(const SortedPoints& second){
			for(size_t i=0,n=this->size();i<n;i++){
				if(Bin(i).X()==second[i].X())
					Bin(i).varY()*=second[i].Y();
				else
					throw Exception<SortedPoints>("Cannot perform arithmetic actions two chains. X coordinates differ");
			}
			return *this;
		}
		SortedPoints&operator*=(const SortedPoints&&second){return operator*=(second);}
		SortedPoints&operator*=(const Func f){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()*=f(Bin(i).X());
			return *this;
		}
		SortedPoints&operator*=(const numY&c){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()*=c;
			return *this;
		}
		SortedPoints&operator*=(const numY&&c){
			return operator*=(c);
		}
		
		SortedPoints&operator/=(const SortedPoints& second){
			for(size_t i=0,n=this->size();i<n;i++){
				if(Bin(i).X()==second[i].X())
					Bin(i).varY()/=second[i].Y();
				else
					throw Exception<SortedPoints>("Cannot perform arithmetic actions two chains. X coordinates differ");
			}
			return *this;
		}
		SortedPoints&operator/=(const SortedPoints&&second){return operator/=(second);}
		SortedPoints&operator/=(const Func f){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()/=f(Bin(i).X());
			return *this;
		}
		SortedPoints&operator/=(const numY&c){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()/=c;
			return *this;
		}
		SortedPoints&operator/=(const numY&&c){
			return operator/=(c);
		}
		
		const SortedPoints operator+(const SortedPoints&other)const{return SortedPoints(*this)+=other;}
		const SortedPoints operator+(const SortedPoints&&other)const{return SortedPoints(*this)+=other;}
		const SortedPoints operator-(const SortedPoints&other)const{return SortedPoints(*this)-=other;}
		const SortedPoints operator-(const SortedPoints&&other)const{return SortedPoints(*this)-=other;}
		const SortedPoints operator*(const SortedPoints&other)const{return SortedPoints(*this)*=other;}
		const SortedPoints operator*(const SortedPoints&&other)const{return SortedPoints(*this)*=other;}
		const SortedPoints operator/(const SortedPoints&other)const{return SortedPoints(*this)/=other;}
		const SortedPoints operator/(const SortedPoints&&other)const{return SortedPoints(*this)/=other;}
		
		const SortedPoints operator+(const numY&other)const{return SortedPoints(*this)+=other;}
		const SortedPoints operator+(const numY&&other)const{return SortedPoints(*this)+=other;}
		const SortedPoints operator-(const numY&other)const{return SortedPoints(*this)-=other;}
		const SortedPoints operator-(const numY&&other)const{return SortedPoints(*this)-=other;}
		const SortedPoints operator*(const numY&other)const{return SortedPoints(*this)*=other;}
		const SortedPoints operator*(const numY&&other)const{return SortedPoints(*this)*=other;}
		const SortedPoints operator/(const numY&other)const{return SortedPoints(*this)/=other;}
		const SortedPoints operator/(const numY&&other)const{return SortedPoints(*this)/=other;}

		const SortedPoints operator+(const Func other)const{return SortedPoints(*this)+=other;}
		const SortedPoints operator-(const Func other)const{return SortedPoints(*this)-=other;}
		const SortedPoints operator*(const Func other)const{return SortedPoints(*this)*=other;}
		const SortedPoints operator/(const Func other)const{return SortedPoints(*this)/=other;}
		
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class BiSortedPoints{
	private:
		SortedChain<numtX> m_x_axis;
		SortedChain<numtY> m_y_axis;
		std::vector<std::vector<numtZ>> m_data;
		void init(){
			m_data.clear();
			for(size_t i=0,I=m_x_axis.size();i<I;i++){
				m_data.push_back(std::vector<numtZ>());
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					m_data[m_data.size()-1].push_back(numtZ(0));
			}
		}
	public:
		const SortedChain<numtX>&X()const{return m_x_axis;}
		const SortedChain<numtY>&Y()const{return m_y_axis;}
		BiSortedPoints(const std::initializer_list<numtX>&X,const std::initializer_list<numtY>&Y)
		:m_x_axis(X),m_y_axis(Y){init();}
		BiSortedPoints(const std::initializer_list<numtX>&&X,const std::initializer_list<numtY>&&Y)
		:m_x_axis(X),m_y_axis(Y){init();}
		BiSortedPoints(const SortedChain<numtX>&X,const SortedChain<numtY>&Y)
		:m_x_axis(X),m_y_axis(Y){init();}
		BiSortedPoints(const SortedChain<numtX>&&X,const SortedChain<numtY>&&Y)
		:m_x_axis(X),m_y_axis(Y){init();}
		BiSortedPoints():BiSortedPoints({},{}){}
		BiSortedPoints(const BiSortedPoints&source):m_x_axis(source.X()),m_y_axis(source.Y()){
			for(size_t i=0,I=source.m_data.size();i<I;i++){
				m_data.push_back(std::vector<numtZ>());
				for(const auto&item:source.m_data[i])
					m_data[i].push_back(item);
			}
		}
		virtual ~BiSortedPoints(){}
		typedef typename std::vector<std::vector<numtZ>>::const_iterator const_iterator;
		const_iterator begin()const{return m_data.cbegin();}
		const_iterator cbegin()const{return m_data.cbegin();}
		const_iterator end() const{return m_data.cend();}
		const_iterator cend() const{return m_data.cend();}
		size_t size()const{return m_data.size();}
		const std::vector<numtZ>&operator[](const size_t i)const{
			if(size()<=i)throw Exception<BiSortedPoints>("range check error");
			return m_data[i];
		}
		numtZ&Bin(const size_t i,const size_t j){
			if(size()<=i)throw Exception<BiSortedPoints>("range check error");
			if(m_data[i].size()<=j)throw Exception<BiSortedPoints>("range check error");
			return m_data[i][j];
		}
		const point3d<numtX,numtY,numtZ> operator()(const size_t i,const size_t j)const{
			if(size()<=i)throw Exception<BiSortedPoints>("range check error");
			if(m_y_axis.size()<=j)throw Exception<BiSortedPoints>("range check error");
			return point3d<numtX,numtY,numtZ>(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycle(const std::function<void(const point3d<numtX,numtY,numtZ>&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++){
					point3d<numtX,numtY,numtZ> P(m_x_axis[i],m_y_axis[j],m_data[i][j]);
					f(P);
				}
		}
		void FullCycle(const std::function<void(const numtX&,const numtY&,const numtZ&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycleVar(const std::function<void(const numtX&,const numtY&,numtZ&)>f){
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		
	};
};
#endif