OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23842731) q[0];
sx q[0];
rz(-2.981346) q[0];
sx q[0];
rz(1.2877553) q[0];
rz(-0.84914485) q[1];
sx q[1];
rz(-1.1916817) q[1];
sx q[1];
rz(-2.3828659) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37739524) q[0];
sx q[0];
rz(-1.8599417) q[0];
sx q[0];
rz(1.3290861) q[0];
rz(-0.24234645) q[2];
sx q[2];
rz(-1.0966612) q[2];
sx q[2];
rz(0.15210064) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0352896) q[1];
sx q[1];
rz(-2.311541) q[1];
sx q[1];
rz(-2.0705877) q[1];
x q[2];
rz(2.4843744) q[3];
sx q[3];
rz(-2.2003678) q[3];
sx q[3];
rz(-2.8384125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3005001) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(0.50660261) q[2];
rz(-1.6182342) q[3];
sx q[3];
rz(-2.5169499) q[3];
sx q[3];
rz(-0.055559572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8697206) q[0];
sx q[0];
rz(-2.6756918) q[0];
sx q[0];
rz(0.063657612) q[0];
rz(1.1977389) q[1];
sx q[1];
rz(-1.5143737) q[1];
sx q[1];
rz(2.1228085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2301529) q[0];
sx q[0];
rz(-0.24556118) q[0];
sx q[0];
rz(-0.6397689) q[0];
rz(-pi) q[1];
rz(-0.34060654) q[2];
sx q[2];
rz(-1.7939027) q[2];
sx q[2];
rz(-2.4702213) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9460822) q[1];
sx q[1];
rz(-1.2942593) q[1];
sx q[1];
rz(2.4184188) q[1];
rz(-1.7131931) q[3];
sx q[3];
rz(-2.3835858) q[3];
sx q[3];
rz(2.6482794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4360875) q[2];
sx q[2];
rz(-2.5477396) q[2];
sx q[2];
rz(-2.2994821) q[2];
rz(1.0574794) q[3];
sx q[3];
rz(-2.2130241) q[3];
sx q[3];
rz(1.9717982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975526) q[0];
sx q[0];
rz(-2.121448) q[0];
sx q[0];
rz(-1.1919588) q[0];
rz(2.2167218) q[1];
sx q[1];
rz(-0.7083188) q[1];
sx q[1];
rz(-1.8398197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98604964) q[0];
sx q[0];
rz(-1.0240004) q[0];
sx q[0];
rz(0.56122551) q[0];
rz(-pi) q[1];
rz(-2.9495732) q[2];
sx q[2];
rz(-2.9387967) q[2];
sx q[2];
rz(1.5206199) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30693248) q[1];
sx q[1];
rz(-1.7494769) q[1];
sx q[1];
rz(2.6684922) q[1];
rz(-pi) q[2];
rz(-2.202192) q[3];
sx q[3];
rz(-2.7223848) q[3];
sx q[3];
rz(-1.4155751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3785582) q[2];
sx q[2];
rz(-3.0578461) q[2];
sx q[2];
rz(-3.0540826) q[2];
rz(1.8267501) q[3];
sx q[3];
rz(-1.0564691) q[3];
sx q[3];
rz(2.1480613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73115504) q[0];
sx q[0];
rz(-1.9696099) q[0];
sx q[0];
rz(-0.80899578) q[0];
rz(-2.1058829) q[1];
sx q[1];
rz(-2.6782942) q[1];
sx q[1];
rz(-2.1083924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9161578) q[0];
sx q[0];
rz(-1.5317129) q[0];
sx q[0];
rz(2.2351511) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43605752) q[2];
sx q[2];
rz(-2.0654709) q[2];
sx q[2];
rz(3.0835754) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.08398) q[1];
sx q[1];
rz(-1.3911493) q[1];
sx q[1];
rz(0.47112005) q[1];
rz(0.71113385) q[3];
sx q[3];
rz(-2.4681849) q[3];
sx q[3];
rz(-2.8720958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0741299) q[2];
sx q[2];
rz(-2.9758657) q[2];
sx q[2];
rz(1.5376252) q[2];
rz(-1.7240723) q[3];
sx q[3];
rz(-1.1105024) q[3];
sx q[3];
rz(-2.0541151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23987016) q[0];
sx q[0];
rz(-2.9481695) q[0];
sx q[0];
rz(1.9523917) q[0];
rz(0.72422782) q[1];
sx q[1];
rz(-1.291357) q[1];
sx q[1];
rz(-1.6667746) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3132202) q[0];
sx q[0];
rz(-0.68183866) q[0];
sx q[0];
rz(0.71558715) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9846813) q[2];
sx q[2];
rz(-2.3381066) q[2];
sx q[2];
rz(-0.29565865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9660419) q[1];
sx q[1];
rz(-1.908424) q[1];
sx q[1];
rz(-0.32159827) q[1];
rz(-pi) q[2];
rz(-0.83051564) q[3];
sx q[3];
rz(-0.6336824) q[3];
sx q[3];
rz(-0.049402852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38179794) q[2];
sx q[2];
rz(-2.2319824) q[2];
sx q[2];
rz(0.46662113) q[2];
rz(-1.4383379) q[3];
sx q[3];
rz(-1.8432901) q[3];
sx q[3];
rz(2.2013825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4284215) q[0];
sx q[0];
rz(-2.3155825) q[0];
sx q[0];
rz(-0.010350479) q[0];
rz(-1.5230644) q[1];
sx q[1];
rz(-1.7489988) q[1];
sx q[1];
rz(-1.7656322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15555602) q[0];
sx q[0];
rz(-1.9080093) q[0];
sx q[0];
rz(2.8854779) q[0];
x q[1];
rz(-2.544246) q[2];
sx q[2];
rz(-0.22946363) q[2];
sx q[2];
rz(-0.84537431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37782323) q[1];
sx q[1];
rz(-1.0727912) q[1];
sx q[1];
rz(-1.3643907) q[1];
rz(-pi) q[2];
rz(2.2155432) q[3];
sx q[3];
rz(-0.9315486) q[3];
sx q[3];
rz(-0.13945068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4202262) q[2];
sx q[2];
rz(-2.0927636) q[2];
sx q[2];
rz(2.4143207) q[2];
rz(0.51491245) q[3];
sx q[3];
rz(-1.8102831) q[3];
sx q[3];
rz(1.4179199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10776831) q[0];
sx q[0];
rz(-1.533778) q[0];
sx q[0];
rz(2.7346101) q[0];
rz(0.96626967) q[1];
sx q[1];
rz(-2.2844908) q[1];
sx q[1];
rz(-3.0028717) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7329368) q[0];
sx q[0];
rz(-1.6870903) q[0];
sx q[0];
rz(-0.035712042) q[0];
x q[1];
rz(1.3499603) q[2];
sx q[2];
rz(-1.2996309) q[2];
sx q[2];
rz(1.3214932) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.87200621) q[1];
sx q[1];
rz(-1.2254189) q[1];
sx q[1];
rz(-0.38265444) q[1];
x q[2];
rz(-1.2267031) q[3];
sx q[3];
rz(-1.4315245) q[3];
sx q[3];
rz(2.1849098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2762642) q[2];
sx q[2];
rz(-1.7636969) q[2];
sx q[2];
rz(0.31065568) q[2];
rz(1.8336512) q[3];
sx q[3];
rz(-1.5111978) q[3];
sx q[3];
rz(-2.7697897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0528316) q[0];
sx q[0];
rz(-2.796687) q[0];
sx q[0];
rz(-3.0533691) q[0];
rz(3.131033) q[1];
sx q[1];
rz(-1.6190395) q[1];
sx q[1];
rz(-0.47058502) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39401606) q[0];
sx q[0];
rz(-0.86578275) q[0];
sx q[0];
rz(0.88253077) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80704524) q[2];
sx q[2];
rz(-1.7084439) q[2];
sx q[2];
rz(1.1040083) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.499191) q[1];
sx q[1];
rz(-2.8676053) q[1];
sx q[1];
rz(-1.8887331) q[1];
rz(-pi) q[2];
rz(1.2180738) q[3];
sx q[3];
rz(-1.439487) q[3];
sx q[3];
rz(-2.2985947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9528902) q[2];
sx q[2];
rz(-2.1850093) q[2];
sx q[2];
rz(1.0799705) q[2];
rz(-2.2522669) q[3];
sx q[3];
rz(-1.5998799) q[3];
sx q[3];
rz(1.5842452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66460669) q[0];
sx q[0];
rz(-2.7598858) q[0];
sx q[0];
rz(2.3241924) q[0];
rz(-2.3951702) q[1];
sx q[1];
rz(-2.4669929) q[1];
sx q[1];
rz(2.2019763) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8259807) q[0];
sx q[0];
rz(-1.5712757) q[0];
sx q[0];
rz(1.5536683) q[0];
rz(-pi) q[1];
rz(-0.65931658) q[2];
sx q[2];
rz(-0.43956471) q[2];
sx q[2];
rz(-2.2826113) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6691362) q[1];
sx q[1];
rz(-1.8207153) q[1];
sx q[1];
rz(-0.12052287) q[1];
x q[2];
rz(-1.1094395) q[3];
sx q[3];
rz(-0.45518866) q[3];
sx q[3];
rz(2.2206642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32999906) q[2];
sx q[2];
rz(-0.96572319) q[2];
sx q[2];
rz(0.3453556) q[2];
rz(-1.6164448) q[3];
sx q[3];
rz(-1.7914146) q[3];
sx q[3];
rz(2.7143872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3578607) q[0];
sx q[0];
rz(-1.1564199) q[0];
sx q[0];
rz(0.082948908) q[0];
rz(-2.9792765) q[1];
sx q[1];
rz(-1.4444618) q[1];
sx q[1];
rz(0.69581318) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37980027) q[0];
sx q[0];
rz(-1.5289991) q[0];
sx q[0];
rz(2.8618852) q[0];
x q[1];
rz(-1.4415353) q[2];
sx q[2];
rz(-1.5737572) q[2];
sx q[2];
rz(-1.2128613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73826927) q[1];
sx q[1];
rz(-1.5898855) q[1];
sx q[1];
rz(2.6973069) q[1];
rz(2.4365594) q[3];
sx q[3];
rz(-1.789973) q[3];
sx q[3];
rz(0.33112835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42056981) q[2];
sx q[2];
rz(-2.9837954) q[2];
sx q[2];
rz(1.921462) q[2];
rz(-0.57341352) q[3];
sx q[3];
rz(-1.8107332) q[3];
sx q[3];
rz(-2.8964608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2777916) q[0];
sx q[0];
rz(-1.5024804) q[0];
sx q[0];
rz(-0.13308751) q[0];
rz(0.0028903891) q[1];
sx q[1];
rz(-0.4160226) q[1];
sx q[1];
rz(-0.70152534) q[1];
rz(-1.3221424) q[2];
sx q[2];
rz(-1.4035618) q[2];
sx q[2];
rz(1.0623912) q[2];
rz(1.5118128) q[3];
sx q[3];
rz(-0.77131693) q[3];
sx q[3];
rz(-0.49280096) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
