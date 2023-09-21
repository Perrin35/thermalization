OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17833248) q[0];
sx q[0];
rz(-1.4890716) q[0];
sx q[0];
rz(-0.89515495) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(-1.4562343) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0460912) q[0];
sx q[0];
rz(-1.7016194) q[0];
sx q[0];
rz(-0.019898947) q[0];
rz(-pi) q[1];
rz(2.7675682) q[2];
sx q[2];
rz(-2.0846539) q[2];
sx q[2];
rz(2.6467269) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1617042) q[1];
sx q[1];
rz(-2.4474505) q[1];
sx q[1];
rz(0.23692268) q[1];
rz(0.31984826) q[3];
sx q[3];
rz(-1.5187129) q[3];
sx q[3];
rz(-0.44714123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75498092) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(2.6888729) q[2];
rz(-2.9833941) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(-0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92900705) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(-1.989495) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-2.6706085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5841056) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(1.4529865) q[0];
x q[1];
rz(-1.2330301) q[2];
sx q[2];
rz(-0.62440364) q[2];
sx q[2];
rz(-2.117346) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4091332) q[1];
sx q[1];
rz(-1.0832548) q[1];
sx q[1];
rz(0.061846102) q[1];
x q[2];
rz(-1.6509634) q[3];
sx q[3];
rz(-0.87682322) q[3];
sx q[3];
rz(0.44192867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(-0.2557959) q[2];
rz(-1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26329041) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.6059748) q[1];
sx q[1];
rz(-0.43930611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8205748) q[0];
sx q[0];
rz(-1.4903755) q[0];
sx q[0];
rz(3.1119425) q[0];
rz(2.7810532) q[2];
sx q[2];
rz(-1.5224824) q[2];
sx q[2];
rz(-2.609032) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9860552) q[1];
sx q[1];
rz(-1.7441161) q[1];
sx q[1];
rz(2.1876513) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4337173) q[3];
sx q[3];
rz(-1.0573514) q[3];
sx q[3];
rz(1.0182667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1674041) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.2444929) q[0];
rz(0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(-0.8746075) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65507209) q[0];
sx q[0];
rz(-0.97210303) q[0];
sx q[0];
rz(2.952851) q[0];
x q[1];
rz(0.38509102) q[2];
sx q[2];
rz(-1.706012) q[2];
sx q[2];
rz(1.1831634) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8352656) q[1];
sx q[1];
rz(-1.3439461) q[1];
sx q[1];
rz(-2.7579685) q[1];
rz(-0.10505418) q[3];
sx q[3];
rz(-1.6451562) q[3];
sx q[3];
rz(-2.8021333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(2.5710035) q[0];
rz(-2.5866306) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(0.74329174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63440454) q[0];
sx q[0];
rz(-0.92538639) q[0];
sx q[0];
rz(-0.27642823) q[0];
rz(-pi) q[1];
rz(-2.6402316) q[2];
sx q[2];
rz(-1.6622346) q[2];
sx q[2];
rz(-2.3102592) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35768269) q[1];
sx q[1];
rz(-2.8257113) q[1];
sx q[1];
rz(2.3481579) q[1];
rz(-pi) q[2];
rz(2.8652142) q[3];
sx q[3];
rz(-2.31782) q[3];
sx q[3];
rz(-2.7793022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1935929) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(2.8732079) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140117) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(-1.0553029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6049833) q[0];
sx q[0];
rz(-0.6491001) q[0];
sx q[0];
rz(-1.8761294) q[0];
rz(-pi) q[1];
rz(-2.2248473) q[2];
sx q[2];
rz(-1.9987717) q[2];
sx q[2];
rz(-0.59631729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1046022) q[1];
sx q[1];
rz(-1.9676419) q[1];
sx q[1];
rz(2.2657822) q[1];
rz(2.8860693) q[3];
sx q[3];
rz(-0.75948411) q[3];
sx q[3];
rz(1.7373191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(-2.4957538) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(-1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(0.77254599) q[0];
rz(1.729471) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(-0.57377446) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6360146) q[0];
sx q[0];
rz(-1.6806707) q[0];
sx q[0];
rz(-1.3789603) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3705809) q[2];
sx q[2];
rz(-1.5491312) q[2];
sx q[2];
rz(-1.3862762) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3473914) q[1];
sx q[1];
rz(-0.83354356) q[1];
sx q[1];
rz(0.79574037) q[1];
x q[2];
rz(-2.9912234) q[3];
sx q[3];
rz(-1.8024076) q[3];
sx q[3];
rz(1.7600972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(-0.77073628) q[2];
rz(-2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.27281) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(-1.1707206) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(-1.2957113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8745236) q[0];
sx q[0];
rz(-2.0282312) q[0];
sx q[0];
rz(-1.9382856) q[0];
rz(-0.92408085) q[2];
sx q[2];
rz(-0.43281049) q[2];
sx q[2];
rz(-0.88976394) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2669275) q[1];
sx q[1];
rz(-1.3473359) q[1];
sx q[1];
rz(-2.2063971) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1300163) q[3];
sx q[3];
rz(-1.5960346) q[3];
sx q[3];
rz(-0.73362918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43508139) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(0.56813017) q[2];
rz(-0.98012296) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(-2.4639159) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(2.303404) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8294551) q[0];
sx q[0];
rz(-1.9337618) q[0];
sx q[0];
rz(-1.8118993) q[0];
rz(-1.6423177) q[2];
sx q[2];
rz(-2.3841249) q[2];
sx q[2];
rz(-2.1069991) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9233401) q[1];
sx q[1];
rz(-2.372101) q[1];
sx q[1];
rz(1.8583276) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25837274) q[3];
sx q[3];
rz(-1.6515454) q[3];
sx q[3];
rz(-1.0563869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(1.8748803) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(-0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72120136) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(-0.23751968) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(2.9097897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3964513) q[0];
sx q[0];
rz(-1.475435) q[0];
sx q[0];
rz(1.061529) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78328697) q[2];
sx q[2];
rz(-0.33595339) q[2];
sx q[2];
rz(-0.17620262) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3540436) q[1];
sx q[1];
rz(-1.429306) q[1];
sx q[1];
rz(0.267412) q[1];
x q[2];
rz(2.9276804) q[3];
sx q[3];
rz(-2.247346) q[3];
sx q[3];
rz(3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1101749) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5768455) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(-2.8167562) q[2];
sx q[2];
rz(-1.3424716) q[2];
sx q[2];
rz(-2.9906103) q[2];
rz(1.590329) q[3];
sx q[3];
rz(-1.3828779) q[3];
sx q[3];
rz(-1.8120017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
