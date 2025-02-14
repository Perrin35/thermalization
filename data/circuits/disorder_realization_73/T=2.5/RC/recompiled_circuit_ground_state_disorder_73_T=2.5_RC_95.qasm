OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(3.3942437) q[0];
sx q[0];
rz(10.630339) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(-1.815058) q[1];
sx q[1];
rz(2.3958652) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691851) q[0];
sx q[0];
rz(-0.23032941) q[0];
sx q[0];
rz(1.3052829) q[0];
x q[1];
rz(-0.029736515) q[2];
sx q[2];
rz(-0.22448891) q[2];
sx q[2];
rz(2.9475074) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89167833) q[1];
sx q[1];
rz(-2.7418156) q[1];
sx q[1];
rz(-2.1348025) q[1];
x q[2];
rz(-2.9336003) q[3];
sx q[3];
rz(-2.0361414) q[3];
sx q[3];
rz(0.057010827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2241609) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(-1.4707627) q[2];
rz(1.6338232) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(-2.2182218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.1976396) q[0];
sx q[0];
rz(1.5684599) q[0];
rz(-0.16954999) q[1];
sx q[1];
rz(-3.0270271) q[1];
sx q[1];
rz(-3.0068908) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41551521) q[0];
sx q[0];
rz(-1.9163791) q[0];
sx q[0];
rz(0.4536566) q[0];
x q[1];
rz(3.1293389) q[2];
sx q[2];
rz(-1.6394221) q[2];
sx q[2];
rz(1.6100281) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4410494) q[1];
sx q[1];
rz(-1.8722539) q[1];
sx q[1];
rz(0.45977199) q[1];
rz(-pi) q[2];
rz(-0.080218519) q[3];
sx q[3];
rz(-1.748198) q[3];
sx q[3];
rz(-0.77840786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1255101) q[2];
sx q[2];
rz(-1.6566015) q[2];
sx q[2];
rz(-0.15277319) q[2];
rz(-1.3656535) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(-0.16807817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.633054) q[0];
sx q[0];
rz(-0.80798739) q[0];
sx q[0];
rz(-0.48164865) q[0];
rz(2.956849) q[1];
sx q[1];
rz(-1.3667204) q[1];
sx q[1];
rz(2.1462323) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6299958) q[0];
sx q[0];
rz(-2.7212486) q[0];
sx q[0];
rz(-2.1612801) q[0];
rz(-pi) q[1];
rz(1.6312509) q[2];
sx q[2];
rz(-1.619307) q[2];
sx q[2];
rz(1.2918351) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0936443) q[1];
sx q[1];
rz(-1.2721095) q[1];
sx q[1];
rz(0.68409749) q[1];
rz(-0.030582436) q[3];
sx q[3];
rz(-0.97506071) q[3];
sx q[3];
rz(0.047614656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92622009) q[2];
sx q[2];
rz(-3.0871349) q[2];
sx q[2];
rz(3.0769297) q[2];
rz(2.0926545) q[3];
sx q[3];
rz(-3.1147396) q[3];
sx q[3];
rz(-1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30508405) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(-2.3205561) q[0];
rz(-0.07846421) q[1];
sx q[1];
rz(-1.4469701) q[1];
sx q[1];
rz(0.99748126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067046384) q[0];
sx q[0];
rz(-0.95854488) q[0];
sx q[0];
rz(-1.9590098) q[0];
x q[1];
rz(1.0216265) q[2];
sx q[2];
rz(-3.0722174) q[2];
sx q[2];
rz(2.3424984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32179896) q[1];
sx q[1];
rz(-2.5410278) q[1];
sx q[1];
rz(2.0482045) q[1];
rz(-pi) q[2];
x q[2];
rz(2.663199) q[3];
sx q[3];
rz(-1.6217303) q[3];
sx q[3];
rz(-2.3277605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0923826) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(-1.1537665) q[2];
rz(2.9554101) q[3];
sx q[3];
rz(-0.081705339) q[3];
sx q[3];
rz(-2.8103099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368197) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(-1.1432884) q[0];
rz(-1.1413057) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(-2.6236261) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2835613) q[0];
sx q[0];
rz(-2.1037344) q[0];
sx q[0];
rz(0.74646797) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5858272) q[2];
sx q[2];
rz(-1.5698729) q[2];
sx q[2];
rz(0.24713384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8570366) q[1];
sx q[1];
rz(-1.2412046) q[1];
sx q[1];
rz(0.28917851) q[1];
x q[2];
rz(-2.84359) q[3];
sx q[3];
rz(-1.4299222) q[3];
sx q[3];
rz(0.44671392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8225857) q[2];
sx q[2];
rz(-0.95390445) q[2];
sx q[2];
rz(0.32773584) q[2];
rz(1.94708) q[3];
sx q[3];
rz(-3.0142398) q[3];
sx q[3];
rz(2.2849042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0026534) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(0.56185454) q[0];
rz(-1.6506763) q[1];
sx q[1];
rz(-1.6249388) q[1];
sx q[1];
rz(3.0432826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54360169) q[0];
sx q[0];
rz(-2.7181135) q[0];
sx q[0];
rz(1.4844358) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19616429) q[2];
sx q[2];
rz(-0.00065302424) q[2];
sx q[2];
rz(3.0214423) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5206388) q[1];
sx q[1];
rz(-0.47157447) q[1];
sx q[1];
rz(0.10271272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6410511) q[3];
sx q[3];
rz(-1.9841474) q[3];
sx q[3];
rz(-2.7374637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0564698) q[2];
sx q[2];
rz(-2.9462908) q[2];
sx q[2];
rz(2.0758212) q[2];
rz(-0.39984518) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(-1.8736418) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14168508) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(2.0211925) q[1];
sx q[1];
rz(-0.1381865) q[1];
sx q[1];
rz(0.33946005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2223245) q[0];
sx q[0];
rz(-1.647729) q[0];
sx q[0];
rz(1.5841106) q[0];
rz(2.5869114) q[2];
sx q[2];
rz(-1.5590057) q[2];
sx q[2];
rz(-1.5761216) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.554018) q[1];
sx q[1];
rz(-1.6563935) q[1];
sx q[1];
rz(-0.0016335131) q[1];
rz(-pi) q[2];
rz(0.33879316) q[3];
sx q[3];
rz(-1.4767891) q[3];
sx q[3];
rz(0.28123048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3642984) q[2];
sx q[2];
rz(-3.1396781) q[2];
sx q[2];
rz(2.778229) q[2];
rz(-1.0904788) q[3];
sx q[3];
rz(-0.57991475) q[3];
sx q[3];
rz(1.1183848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5217487) q[0];
sx q[0];
rz(-2.813297) q[0];
sx q[0];
rz(-0.92292619) q[0];
rz(-1.482831) q[1];
sx q[1];
rz(-0.62983477) q[1];
sx q[1];
rz(0.0042075687) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370458) q[0];
sx q[0];
rz(-2.1677289) q[0];
sx q[0];
rz(-2.202522) q[0];
x q[1];
rz(3.1356631) q[2];
sx q[2];
rz(-1.9816795) q[2];
sx q[2];
rz(3.1290999) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6018659) q[1];
sx q[1];
rz(-1.5063573) q[1];
sx q[1];
rz(2.0468725) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0757853) q[3];
sx q[3];
rz(-0.61344693) q[3];
sx q[3];
rz(-0.25019161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5620455) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(1.2008249) q[2];
rz(-1.3935401) q[3];
sx q[3];
rz(-0.0037007185) q[3];
sx q[3];
rz(-2.4414731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30018184) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(1.9218943) q[0];
rz(1.3297184) q[1];
sx q[1];
rz(-2.0025573) q[1];
sx q[1];
rz(-2.9776998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4224098) q[0];
sx q[0];
rz(-1.2162195) q[0];
sx q[0];
rz(2.9585881) q[0];
x q[1];
rz(0.58456011) q[2];
sx q[2];
rz(-1.341408) q[2];
sx q[2];
rz(3.0363415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6299705) q[1];
sx q[1];
rz(-0.13221201) q[1];
sx q[1];
rz(2.0396784) q[1];
x q[2];
rz(-2.2396062) q[3];
sx q[3];
rz(-0.94873442) q[3];
sx q[3];
rz(-0.40611503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4280052) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(-3.0981787) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49122214) q[0];
sx q[0];
rz(-0.0064042052) q[0];
sx q[0];
rz(0.48625913) q[0];
rz(-0.69475118) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(-0.4756701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.75695) q[0];
sx q[0];
rz(-1.5371377) q[0];
sx q[0];
rz(1.606694) q[0];
x q[1];
rz(2.2769502) q[2];
sx q[2];
rz(-1.5967895) q[2];
sx q[2];
rz(3.1033422) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.039363843) q[1];
sx q[1];
rz(-0.96323955) q[1];
sx q[1];
rz(1.5629014) q[1];
rz(-0.31020152) q[3];
sx q[3];
rz(-1.3942914) q[3];
sx q[3];
rz(1.4909397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67705578) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(-1.2976868) q[2];
rz(-1.248598) q[3];
sx q[3];
rz(-2.5864351) q[3];
sx q[3];
rz(-2.7738074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0155335) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(0.88232782) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(-1.6068983) q[2];
sx q[2];
rz(-1.6629728) q[2];
sx q[2];
rz(-2.865553) q[2];
rz(-0.020261665) q[3];
sx q[3];
rz(-0.23923418) q[3];
sx q[3];
rz(-3.1399865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
