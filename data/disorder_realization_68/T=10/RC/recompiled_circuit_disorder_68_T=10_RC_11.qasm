OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7473937) q[0];
sx q[0];
rz(-2.6497901) q[0];
sx q[0];
rz(2.9536182) q[0];
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(2.7741073) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1023094) q[0];
sx q[0];
rz(-1.6741721) q[0];
sx q[0];
rz(-2.1084059) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0796702) q[2];
sx q[2];
rz(-1.6899781) q[2];
sx q[2];
rz(-1.0246547) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5012813) q[1];
sx q[1];
rz(-0.2959364) q[1];
sx q[1];
rz(-2.179115) q[1];
rz(-1.7484619) q[3];
sx q[3];
rz(-1.2347617) q[3];
sx q[3];
rz(0.48776585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.964103) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(2.5906079) q[2];
rz(1.8356813) q[3];
sx q[3];
rz(-1.4923613) q[3];
sx q[3];
rz(-1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47857639) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(-0.42981237) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(-2.205251) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443003) q[0];
sx q[0];
rz(-2.8951277) q[0];
sx q[0];
rz(-0.36578567) q[0];
rz(-1.9905375) q[2];
sx q[2];
rz(-2.310576) q[2];
sx q[2];
rz(-1.8120399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.7649819) q[1];
sx q[1];
rz(-0.88135834) q[1];
sx q[1];
rz(-1.7060243) q[1];
rz(-pi) q[2];
rz(2.0735047) q[3];
sx q[3];
rz(-1.5912676) q[3];
sx q[3];
rz(-2.3308144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77461809) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(0.42638391) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(-3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8957829) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(0.93908969) q[0];
rz(2.242873) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(-0.59392196) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1840738) q[0];
sx q[0];
rz(-0.44555095) q[0];
sx q[0];
rz(-0.81934388) q[0];
x q[1];
rz(1.5005641) q[2];
sx q[2];
rz(-2.4519081) q[2];
sx q[2];
rz(0.94674142) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18332874) q[1];
sx q[1];
rz(-1.7839583) q[1];
sx q[1];
rz(-2.0002685) q[1];
rz(-pi) q[2];
rz(-2.5251758) q[3];
sx q[3];
rz(-2.5865002) q[3];
sx q[3];
rz(-2.4544231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64017355) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(1.7017986) q[2];
rz(0.38763186) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(2.7534527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(0.50278062) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(0.75685135) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8096136) q[0];
sx q[0];
rz(-2.6410714) q[0];
sx q[0];
rz(2.39397) q[0];
x q[1];
rz(-1.5449764) q[2];
sx q[2];
rz(-0.46586793) q[2];
sx q[2];
rz(2.4398838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13285747) q[1];
sx q[1];
rz(-0.99616226) q[1];
sx q[1];
rz(0.43032129) q[1];
rz(0.54550708) q[3];
sx q[3];
rz(-0.81199284) q[3];
sx q[3];
rz(-3.0363887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7148774) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(1.654401) q[2];
rz(2.5590844) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(-2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(-1.2879397) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(-2.0910738) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30848635) q[0];
sx q[0];
rz(-1.3163438) q[0];
sx q[0];
rz(1.2480877) q[0];
rz(0.96111416) q[2];
sx q[2];
rz(-0.69176199) q[2];
sx q[2];
rz(1.320653) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.037462385) q[1];
sx q[1];
rz(-2.3298652) q[1];
sx q[1];
rz(-0.8123668) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53253048) q[3];
sx q[3];
rz(-0.88056394) q[3];
sx q[3];
rz(1.1158451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.918255) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(-2.5081432) q[2];
rz(-1.9472306) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(2.8884086) q[0];
rz(-1.5340012) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(1.6794499) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40186858) q[0];
sx q[0];
rz(-1.3603633) q[0];
sx q[0];
rz(-1.6519288) q[0];
x q[1];
rz(0.92163779) q[2];
sx q[2];
rz(-1.3772794) q[2];
sx q[2];
rz(-0.18093872) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6219382) q[1];
sx q[1];
rz(-2.7192273) q[1];
sx q[1];
rz(2.7154891) q[1];
x q[2];
rz(1.2268279) q[3];
sx q[3];
rz(-1.868639) q[3];
sx q[3];
rz(0.96836585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9399461) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(-0.80491006) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-1.0837519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2729623) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(-1.0246798) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(2.129508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597848) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(-2.3083789) q[0];
rz(-1.4321248) q[2];
sx q[2];
rz(-1.9562634) q[2];
sx q[2];
rz(-3.1388381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.183179) q[1];
sx q[1];
rz(-1.603754) q[1];
sx q[1];
rz(-2.5977913) q[1];
rz(-pi) q[2];
rz(0.13109644) q[3];
sx q[3];
rz(-2.5790865) q[3];
sx q[3];
rz(-1.2870115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6107789) q[2];
sx q[2];
rz(-1.4799708) q[2];
sx q[2];
rz(-2.7116595) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(-2.608192) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124509) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(-2.9274143) q[0];
rz(-2.0902436) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(2.8578551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8003214) q[0];
sx q[0];
rz(-0.30297908) q[0];
sx q[0];
rz(0.11462258) q[0];
rz(-pi) q[1];
rz(-1.0546513) q[2];
sx q[2];
rz(-2.7661341) q[2];
sx q[2];
rz(1.5309389) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.022790837) q[1];
sx q[1];
rz(-2.6036501) q[1];
sx q[1];
rz(-2.7522037) q[1];
rz(2.22105) q[3];
sx q[3];
rz(-1.0245819) q[3];
sx q[3];
rz(-1.2342681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45067898) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(1.696375) q[2];
rz(-1.5971659) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63672367) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(0.069256393) q[0];
rz(-1.6537369) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(-1.5725296) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42800316) q[0];
sx q[0];
rz(-1.1134976) q[0];
sx q[0];
rz(-2.2767115) q[0];
x q[1];
rz(-0.32585085) q[2];
sx q[2];
rz(-1.648765) q[2];
sx q[2];
rz(0.98380145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3410586) q[1];
sx q[1];
rz(-1.3375999) q[1];
sx q[1];
rz(-1.6569767) q[1];
x q[2];
rz(-2.9326354) q[3];
sx q[3];
rz(-0.93271241) q[3];
sx q[3];
rz(-2.0861422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1853603) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(0.64176732) q[0];
rz(1.2311252) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(2.8737601) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5726686) q[0];
sx q[0];
rz(-0.84416443) q[0];
sx q[0];
rz(0.6092682) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1688813) q[2];
sx q[2];
rz(-2.9209666) q[2];
sx q[2];
rz(2.0711183) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47083449) q[1];
sx q[1];
rz(-1.8948312) q[1];
sx q[1];
rz(1.2653989) q[1];
rz(0.28585163) q[3];
sx q[3];
rz(-1.0470069) q[3];
sx q[3];
rz(1.3956192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81007593) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(-0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(0.77990445) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(0.71074625) q[2];
sx q[2];
rz(-1.1087316) q[2];
sx q[2];
rz(1.0089594) q[2];
rz(-1.7846617) q[3];
sx q[3];
rz(-0.90448096) q[3];
sx q[3];
rz(2.7703551) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];