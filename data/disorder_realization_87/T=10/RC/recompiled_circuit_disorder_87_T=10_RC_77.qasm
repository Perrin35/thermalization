OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(2.6013241) q[1];
sx q[1];
rz(7.2202914) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9286081) q[0];
sx q[0];
rz(-2.1323418) q[0];
sx q[0];
rz(-1.0667849) q[0];
rz(-2.5295528) q[2];
sx q[2];
rz(-1.1571552) q[2];
sx q[2];
rz(0.94758247) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.571849) q[1];
sx q[1];
rz(-1.3880001) q[1];
sx q[1];
rz(3.0621431) q[1];
x q[2];
rz(0.097221656) q[3];
sx q[3];
rz(-1.8386278) q[3];
sx q[3];
rz(0.83084805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(3.0554331) q[2];
rz(2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(2.0479726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57698292) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.3264867) q[0];
rz(1.2558698) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(2.870141) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4685681) q[0];
sx q[0];
rz(-0.66267555) q[0];
sx q[0];
rz(2.2492692) q[0];
x q[1];
rz(-2.62918) q[2];
sx q[2];
rz(-2.5645442) q[2];
sx q[2];
rz(-1.3702099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27122341) q[1];
sx q[1];
rz(-1.1961812) q[1];
sx q[1];
rz(0.80925525) q[1];
rz(-pi) q[2];
rz(0.80941697) q[3];
sx q[3];
rz(-1.8339694) q[3];
sx q[3];
rz(2.0078878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0469971) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(0.57717741) q[2];
rz(-2.2180637) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975824) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(2.999021) q[0];
rz(1.3525195) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(2.9325063) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81461834) q[0];
sx q[0];
rz(-1.0083535) q[0];
sx q[0];
rz(2.4787865) q[0];
rz(0.44720165) q[2];
sx q[2];
rz(-2.3502091) q[2];
sx q[2];
rz(0.42391047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28029728) q[1];
sx q[1];
rz(-2.0164844) q[1];
sx q[1];
rz(-2.6764826) q[1];
rz(-pi) q[2];
x q[2];
rz(2.530982) q[3];
sx q[3];
rz(-2.7770677) q[3];
sx q[3];
rz(-0.025346905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7769988) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(-0.40412942) q[2];
rz(-1.2858307) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(0.16734853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.4362815) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(0.96281111) q[0];
rz(0.46936938) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2179746) q[0];
sx q[0];
rz(-1.0110564) q[0];
sx q[0];
rz(0.864242) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5702815) q[2];
sx q[2];
rz(-1.4467903) q[2];
sx q[2];
rz(0.73881432) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.094593781) q[1];
sx q[1];
rz(-1.7013036) q[1];
sx q[1];
rz(2.1427437) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13929316) q[3];
sx q[3];
rz(-0.5738429) q[3];
sx q[3];
rz(-1.3822671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0756388) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(-3.0299419) q[2];
rz(-2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951293) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(2.7850889) q[0];
rz(2.6351392) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(-2.8809663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.168419) q[0];
sx q[0];
rz(-1.6126313) q[0];
sx q[0];
rz(0.14194685) q[0];
x q[1];
rz(1.3172651) q[2];
sx q[2];
rz(-1.8142482) q[2];
sx q[2];
rz(1.4075116) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.063559859) q[1];
sx q[1];
rz(-1.5531925) q[1];
sx q[1];
rz(-0.5743919) q[1];
rz(1.4277677) q[3];
sx q[3];
rz(-1.4266532) q[3];
sx q[3];
rz(-0.10728697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9115209) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(-0.22932209) q[2];
rz(-2.5991332) q[3];
sx q[3];
rz(-2.8312603) q[3];
sx q[3];
rz(1.0361766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(1.0391327) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6468069) q[0];
sx q[0];
rz(-1.8456869) q[0];
sx q[0];
rz(1.5216212) q[0];
rz(-pi) q[1];
rz(0.44712375) q[2];
sx q[2];
rz(-0.56338718) q[2];
sx q[2];
rz(-2.2396357) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9767178) q[1];
sx q[1];
rz(-0.70189447) q[1];
sx q[1];
rz(-0.81275069) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7870951) q[3];
sx q[3];
rz(-1.6445451) q[3];
sx q[3];
rz(-2.4019965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.002939) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(0.49267832) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7572927) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(-3.0199155) q[0];
rz(-1.9901468) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(2.7391403) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84797317) q[0];
sx q[0];
rz(-1.3012039) q[0];
sx q[0];
rz(-0.86975354) q[0];
x q[1];
rz(1.5624814) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(2.7260821) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35343364) q[1];
sx q[1];
rz(-2.806059) q[1];
sx q[1];
rz(-0.56281705) q[1];
rz(-pi) q[2];
rz(1.6472858) q[3];
sx q[3];
rz(-0.34919958) q[3];
sx q[3];
rz(-2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1272614) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(0.86223117) q[2];
rz(-2.6640653) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(-2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(-0.73295897) q[0];
rz(0.14006242) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(1.0345116) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1326133) q[0];
sx q[0];
rz(-1.7404557) q[0];
sx q[0];
rz(-1.3404113) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60136749) q[2];
sx q[2];
rz(-0.12430087) q[2];
sx q[2];
rz(2.3483495) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4701177) q[1];
sx q[1];
rz(-0.56709328) q[1];
sx q[1];
rz(1.6398029) q[1];
rz(1.6789867) q[3];
sx q[3];
rz(-1.9868317) q[3];
sx q[3];
rz(-3.0707404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(0.77587664) q[2];
rz(2.4173229) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(-1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6999321) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(0.016816703) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(0.7787849) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6803857) q[0];
sx q[0];
rz(-1.9068423) q[0];
sx q[0];
rz(1.3075605) q[0];
rz(-pi) q[1];
rz(1.0672827) q[2];
sx q[2];
rz(-2.6138966) q[2];
sx q[2];
rz(2.6380981) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89214954) q[1];
sx q[1];
rz(-1.3378694) q[1];
sx q[1];
rz(-1.1942785) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.064568297) q[3];
sx q[3];
rz(-1.7913006) q[3];
sx q[3];
rz(1.8622423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-2.4251535) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.107782) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-0.15144908) q[0];
rz(-0.17627136) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(2.418628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6930981) q[0];
sx q[0];
rz(-1.7211595) q[0];
sx q[0];
rz(3.0340414) q[0];
rz(-1.1935913) q[2];
sx q[2];
rz(-1.7282439) q[2];
sx q[2];
rz(2.8709656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1180229) q[1];
sx q[1];
rz(-1.5155063) q[1];
sx q[1];
rz(-1.0667332) q[1];
x q[2];
rz(2.74182) q[3];
sx q[3];
rz(-2.3120566) q[3];
sx q[3];
rz(-1.5012036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4621949) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(-0.52552137) q[2];
rz(0.28371352) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(0.39653683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59239607) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(-1.0992959) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(2.1500258) q[2];
sx q[2];
rz(-0.75800037) q[2];
sx q[2];
rz(-2.6108685) q[2];
rz(0.59393926) q[3];
sx q[3];
rz(-2.754302) q[3];
sx q[3];
rz(2.1828628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
