//Parameters:
//int polymer: the polymer chosen to be generated
//int position: the cite that the new bond is attached to
void BondInitial(int polymer,int position)
{
int bondx,bondy; // coordinate of bond vectors
int directx,directy;  //directions of bond vectors
int bondlength = 14;
// generate a bond restricted to the set B
while ((bondlength > 13)||(bondlength < 4) )
{
// randomly choose lengths
bondx = (int)(4*RandomFunc2(&Seed));
bondy = (int)(4*RandomFunc2(&Seed));
bondlength = bondx*bondx + bondy*bondy;
if ((bondlength <= 13) &&(bondlength >= 4))
{
// randomly choose directions
if (RandomFunc2(&Seed)>0.5)
directx =-1;
else
directx =1;
if (RandomFunc2(&Seed)>0.5)
directy =-1;
else
directy =1;
// add the new bond to the existing segments
Monomer[polymer][position].MonomerX
= Monomer[polymer][position-1].MonomerX + directx*bondx;
Monomer[polymer][position].MonomerY
= Monomer[polymer][position-1].MonomerY + directy*bondy;
// store the information of the new bond
Monomer[polymer][position].MonomerFormerX
= Monomer[polymer][position].MonomerX;
Monomer[polymer][position].MonomerFormerY
= Monomer[polymer][position].MonomerY;
} }
}

int BondInitialChecking(int polymer, int position)
{
int     BondLength;
double  NewBondMidX, NewBondMidY, NewBondMidZ, OldBondMidX,OldBondMidY
        , OldBondMidZ;
double  CheckLengthX, CheckLengthY, CheckLengthZ, CheckLength;
// check self avoidance
for(int p=0; p < position;p++)
{
// calculate the square of the monomer-monomer distance
 BondLength =
 (Monomer[polymer][position].MonomerX - Monomer[polymer][p].MonomerX)
 *(Monomer[polymer][position].MonomerX-Monomer[polymer][p].MonomerX)
+(Monomer[polymer][position].MonomerY - Monomer[polymer][p].MonomerY)
*(Monomer[polymer][position].MonomerY-Monomer[polymer][p].MonomerY);
// calculate the square of the center-to-center distance
 NewBondMidX= (Monomer[polymer][position].MonomerX +
              Monomer[polymer][position-1].MonomerX)/2.0;
 NewBondMidY= (Monomer[polymer][position].MonomerY +
              Monomer[polymer][position-1].MonomerY)/2.0;

if(p!=0)
{
OldBondMidX =
(Monomer[polymer][p].MonomerX + Monomer[polymer][p-1].MonomerX)/2.0;
OldBondMidY =
(Monomer[polymer][p].MonomerY + Monomer[polymer][p-1].MonomerY)/2.0;
CheckLengthX= NewBondMidX - OldBondMidX;
CheckLengthY= NewBondMidY - OldBondMidY;
CheckLength = CheckLengthX*CheckLengthX + CheckLengthY*CheckLengthY;
}
 else
CheckLength=4;
// check the self-avoiding conditions
if((BondLength < 4)||(CheckLength<2))
 return(0);
 }
 return(1);
}

Parameters:
int polymer: the index of polymer
int position: the monomer that is chosen to move
int choice:  = 1, move x coordinate; =2, move y coordinate
int direction: the direction of the move ( = 1 or -1)
void MonomerMoving(int polymer,int position,int choice,int direction)
{
// keep the information of the original chain for fear that the new
   move disobeys the self avoidance
Monomer[polymer][position].MonomerFormerX
= Monomer[polymer][position].MonomerX;
Monomer[polymer][position].MonomerFormerY
= Monomer[polymer][position].MonomerY;
// choose which coordinate to move
 switch(choice)
{
// moving x coordinate
case 1:
Monomer[polymer][position].MonomerX
= Monomer[polymer][position].MonomerX + direction;
Monomer[polymer][position].MonomerY
= Monomer[polymer][position].MonomerY;
break;
// moving y coordinate
case 2:
Monomer[polymer][position].MonomerX
= Monomer[polymer][position].MonomerX;
Monomer[polymer][position].MonomerY
= Monomer[polymer][position].MonomerY + direction;
break;
default:
break; }
}

Parameters:
int polymer: the index of polymers
int site: the monomer that needs to check
int ran: a random number
int JudgeNewBond(int polymer,int site,double ran )
{
int BondLengthSqu=0,NeighbourSite=0,NeighbourPolymer=0,Temp=0;
int Coordinate1X=0, Coordinate1Y=0, Coordinate2X=0, Coordinate2Y=0;
int Coordinate3X=0, Coordinate3Y=0, Coordinate0X=0, Coordinate0Y=0;
double BendingFirstEnergy=0,BendingSecondEnergy=0,BendingThirdEnergy=0
       , ForceEnergyChange=0;
double randomnumber;
int A0=0,A1=0,A2=0,A3=0;
ForceEnergyChange=0;
randomnumber=ran;
Temp=polymer*NumMonomer+site;
// initialize the bending energy; for flexible polymer, it equals to 0
Monomer[polymer][site].BendingFormerEnergy
=Monomer[polymer][site].BendingEnergy;
Monomer[polymer][site].BondLengthFormerSquare
=Monomer[polymer][site].BondLengthSquare;
//  the following calculates the bending energy, and accept the new
    configuration with the Metropolis acceptance rate
if (site==1)
{
Coordinate3X
= Monomer[polymer][site+2].MonomerX-Monomer[polymer][site+1].MonomerX;
Coordinate3Y
= Monomer[polymer][site+2].MonomerY-Monomer[polymer][site+1].MonomerY;
A3
= Coordinate3X*Coordinate3X+Coordinate3Y*Coordinate3Y;
BendingSecondEnergy = Stiffness
*(1-(Coordinate2X*Coordinate1X+Coordinate2Y*Coordinate1Y)/sqrt(A2*A1));
BendingThirdEnergy =  Stiffness
*(1-(Coordinate2X*Coordinate3X+Coordinate2Y*Coordinate3Y)/sqrt(A2*A3));
// Metropolis acceptance rate
if(exp(BendingSecondEnergy+BendingThirdEnergy
-Monomer[polymer][site+1].BendingEnergy
-Monomer[polymer][site].BendingEnergy+ForceEnergyChange)<randomnumber)
return(0);
 // new configuration fails in the check; recovers the old information
Monomer[polymer][site].BendingEnergy=BendingSecondEnergy;
Monomer[polymer][site+1].BendingEnergy=BendingThirdEnergy;
Monomer[polymer][site].BondLengthSquare=A1;
Monomer[polymer][site+1].BondLengthSquare=A2;
}
// same procedure for site=(NumMonomer-2)
else if(site==(NumMonomer-2))
{
Coordinate0X
= Monomer[polymer][site-1].MonomerX-Monomer[polymer][site-2].MonomerX;
Coordinate0Y
= Monomer[polymer][site-1].MonomerY-Monomer[polymer][site-2].MonomerY;
A0 = Coordinate0X*Coordinate0X+Coordinate0Y*Coordinate0Y;
BendingSecondEnergy =Stiffness*
(1-(Coordinate2X*Coordinate1X+Coordinate2Y*Coordinate1Y)/sqrt(A2*A1));
BendingFirstEnergy =  Stiffness*
(1-(Coordinate1X*Coordinate0X+Coordinate1Y*Coordinate0Y)/sqrt(A1*A0));
if(exp(BendingSecondEnergy+BendingFirstEnergy
-Monomer[polymer][site].BendingEnergy
-Monomer[polymer][site-1].BendingEnergy+ForceEnergyChange)<randomnumber)
return(0);
Monomer[polymer][site].BendingEnergy=BendingSecondEnergy;
Monomer[polymer][site-1].BendingEnergy=BendingFirstEnergy;
Monomer[polymer][site].BondLengthSquare=A1;
Monomer[polymer][site+1].BondLengthSquare=A2;
 }
// same procedure for other sites
else if((site!=(NumMonomer-2))&&(site!=1))
{
Coordinate0X = Monomer[polymer][site-1].MonomerX-Monomer[polymer][site-2].MonomerX;
Coordinate0Y = Monomer[polymer][site-1].MonomerY-Monomer[polymer][site-2].MonomerY;
Coordinate3X = Monomer[polymer][site+2].MonomerX-Monomer[polymer][site+1].MonomerX;
Coordinate3Y = Monomer[polymer][site+2].MonomerY-Monomer[polymer][site+1].MonomerY;
A0 = Coordinate0X*Coordinate0X+Coordinate0Y*Coordinate0Y;
A3 = Coordinate3X*Coordinate3X+Coordinate3Y*Coordinate3Y;
BendingSecondEnergy=Stiffness*
(1-(Coordinate2X*Coordinate1X +Coordinate2Y*Coordinate1Y)/sqrt(A2*A1));
BendingFirstEnergy =  Stiffness*
(1-(Coordinate1X*Coordinate0X+Coordinate1Y*Coordinate0Y)/sqrt(A1*A0));
BendingThirdEnergy =  Stiffness*
(1-(Coordinate2X*Coordinate3X+Coordinate2Y*Coordinate3Y)/sqrt(A2*A3));
if(exp(BendingSecondEnergy+BendingFirstEnergy+BendingThirdEnergy
-Monomer[polymer][site+1].BendingEnergy
-Monomer[polymer][site].BendingEnergy
-Monomer[polymer][site-1].BendingEnergy
+ForceEnergyChange)<randomnumber)
return(0);
Monomer[polymer][site].BendingEnergy=BendingSecondEnergy;
Monomer[polymer][site-1].BendingEnergy=BendingFirstEnergy;
Monomer[polymer][site+1].BendingEnergy=BendingThirdEnergy;
Monomer[polymer][site].BondLengthSquare=A1;
Monomer[polymer][site+1].BondLengthSquare=A2;

}
// calculate the new bond length for site = 0
if(site==0)
{
Coordinate2X
= Monomer[polymer][site+1].MonomerX - Monomer[polymer][site].MonomerX;
Coordinate2Y
= Monomer[polymer][site+1].MonomerY - Monomer[polymer][site].MonomerY;
Coordinate3X
= Monomer[polymer][site+2].MonomerX - Monomer[polymer][site+1].MonomerX;
Coordinate3Y
= Monomer[polymer][site+2].MonomerY - Monomer[polymer][site+1].MonomerY;
A2= Coordinate2X*Coordinate2X + Coordinate2Y*Coordinate2Y;
A3= Coordinate3X*Coordinate3X + Coordinate3Y*Coordinate3Y;
// break the self avoidance
if ((A2<4)||(A2>13))
return(0);
//calculate the monomer-monomer distance between different chains
for(int p=0;p < NumPolymer;p++)
{
if(p == polymer)
{
for(int m=2;m<NumMonomer;m++)
{
BondLengthSqu
= (Monomer[polymer][site].MonomerX - Monomer[p][m].MonomerX)
  *(Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
  + (Monomer[polymer][site].MonomerY - Monomer[p][m].MonomerY)
  *(Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY);
if (BondLengthSqu<4)
return(0);
}
} else
{ for(int m=0;m<NumMonomer;m++)
{
BondLengthSqu
= (Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
   *(Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
   +(Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY)
   *(Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY);
 if (BondLengthSqu<4)
return(0);
} }
}

return(1); }
// the same procedure as site=0; the calculations are slightly different
  for the different position of new monomer
else if (site==(NumMonomer-1))
{
Coordinate1X
= Monomer[polymer][site].MonomerX-Monomer[polymer][site-1].MonomerX;
Coordinate1Y
= Monomer[polymer][site].MonomerY-Monomer[polymer][site-1].MonomerY;
Coordinate0X
= Monomer[polymer][site-1].MonomerX-Monomer[polymer][site-2].MonomerX;
Coordinate0Y
= Monomer[polymer][site-1].MonomerY-Monomer[polymer][site-2].MonomerY;
A0= Coordinate0X*Coordinate0X + Coordinate0Y*Coordinate0Y;
A1= Coordinate1X*Coordinate1X + Coordinate1Y*Coordinate1Y;
if ((A1 < 4)||(A1 > 13))
return(0);
 for(int p=0;p<NumPolymer;p++)
{
if(p == polymer)
{
for(int m=0;m<site;m++)
{
BondLengthSqu
= (Monomer[polymer][site].MonomerX - Monomer[p][m].MonomerX)
  *(Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
  +(Monomer[polymer][site].MonomerY - Monomer[p][m].MonomerY)
   *(Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY);
if (BondLengthSqu<4)
return(0);
}
}
else
{
for(int m=0;m<NumMonomer;m++)
{
BondLengthSqu=
(Monomer[polymer][site].MonomerX - Monomer[p][m].MonomerX)
*(Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
+(Monomer[polymer][site].MonomerY - Monomer[p][m].MonomerY)
*(Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY);
 if (BondLengthSqu<4)
 return(0);
 }
} }
Monomer[polymer][site-1].BendingEnergy=BendingFirstEnergy;

Monomer[polymer][site].BondLengthSquare=A0;
}
else if ((site!=(NumMonomer-1))&&(site!=0))
{
Coordinate1X
= Monomer[polymer][site].MonomerX - Monomer[polymer][site-1].MonomerX;
Coordinate1Y
= Monomer[polymer][site].MonomerY - Monomer[polymer][site-1].MonomerY;
Coordinate2X
= Monomer[polymer][site+1].MonomerX - Monomer[polymer][site].MonomerX;
Coordinate2Y
= Monomer[polymer][site+1].MonomerY - Monomer[polymer][site].MonomerY;
A1= Coordinate1X*Coordinate1X + Coordinate1Y*Coordinate1Y;
A2= Coordinate2X*Coordinate2X + Coordinate2Y*Coordinate2Y;
if ((A1<4)||(A1>13))
return(0);
if ((A2<4)||(A2>13))
return(0);
for(int p=0;p<NumPolymer;p++)
{
if(p==polymer)
{
for(int m=0;m<NumMonomer;m++)
{
if(m!=site)
{
BondLengthSqu
= (Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
  *(Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
  + (Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY)
  *(Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY);
if (BondLengthSqu<4)
return(0);
}
}
}
else
{
for(int m=0;m<NumMonomer;m++)
BondLengthSqu
= (Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
*(Monomer[polymer][site].MonomerX-Monomer[p][m].MonomerX)
+ (Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY)
*(Monomer[polymer][site].MonomerY-Monomer[p][m].MonomerY);
if (BondLengthSqu<4)
return(0);
} }
} }
return(1); }

Parameters:
int order: =0 for the case of two chains during the initialization;
           otherwise =1
int polymer: index of polymers
int X, int Y: x and y coordinate of monomer
int JudgeBoundray(int order,int polymer,int X,int Y)
{
// x and y coordinate are limited in the confinement
 if((X < 0)||(X > (int)(Box_size_x*2.8)))
  return(0);
 if((Y < 0)||(Y> (int)(Box_size_y*2.8)))
     return(0);
// the case of two chains during initialization
 if(order == 0)
  {
    if((Y <= (int)(2.8*Box_size_y/NumPolymer*polymer))
||(Y >= (int)(2.8*Box_size_y/NumPolymer*(polymer+1))))
return(0); }
return(1); }

