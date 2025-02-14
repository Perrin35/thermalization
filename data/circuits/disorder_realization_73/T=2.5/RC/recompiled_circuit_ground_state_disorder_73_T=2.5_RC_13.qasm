OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8483228) q[0];
sx q[0];
rz(-0.25265101) q[0];
sx q[0];
rz(-1.2055612) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(4.4681273) q[1];
sx q[1];
rz(11.820643) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4843895) q[0];
sx q[0];
rz(-1.5108539) q[0];
sx q[0];
rz(-1.3482698) q[0];
rz(1.5640075) q[2];
sx q[2];
rz(-1.3464084) q[2];
sx q[2];
rz(0.16358384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2499143) q[1];
sx q[1];
rz(-0.39977705) q[1];
sx q[1];
rz(2.1348025) q[1];
rz(1.9609591) q[3];
sx q[3];
rz(-0.50658617) q[3];
sx q[3];
rz(0.49662874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2241609) q[2];
sx q[2];
rz(-1.5459583) q[2];
sx q[2];
rz(-1.6708299) q[2];
rz(-1.5077695) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(0.9233709) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(1.5731328) q[0];
rz(-2.9720427) q[1];
sx q[1];
rz(-3.0270271) q[1];
sx q[1];
rz(3.0068908) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54808211) q[0];
sx q[0];
rz(-2.5786886) q[0];
sx q[0];
rz(-2.4538732) q[0];
rz(3.1293389) q[2];
sx q[2];
rz(-1.5021706) q[2];
sx q[2];
rz(1.5315646) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0162279) q[1];
sx q[1];
rz(-2.008359) q[1];
sx q[1];
rz(1.9047649) q[1];
rz(-0.080218519) q[3];
sx q[3];
rz(-1.3933946) q[3];
sx q[3];
rz(0.77840786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1255101) q[2];
sx q[2];
rz(-1.6566015) q[2];
sx q[2];
rz(-0.15277319) q[2];
rz(1.3656535) q[3];
sx q[3];
rz(-3.1047265) q[3];
sx q[3];
rz(2.9735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50853866) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(2.659944) q[0];
rz(0.18474361) q[1];
sx q[1];
rz(-1.3667204) q[1];
sx q[1];
rz(-2.1462323) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99670519) q[0];
sx q[0];
rz(-1.9166244) q[0];
sx q[0];
rz(-0.24390999) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89389385) q[2];
sx q[2];
rz(-0.077493103) q[2];
sx q[2];
rz(2.187196) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75835704) q[1];
sx q[1];
rz(-0.92229533) q[1];
sx q[1];
rz(1.1926257) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.030582436) q[3];
sx q[3];
rz(-0.97506071) q[3];
sx q[3];
rz(0.047614656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2153726) q[2];
sx q[2];
rz(-3.0871349) q[2];
sx q[2];
rz(0.064662956) q[2];
rz(-1.0489382) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30508405) q[0];
sx q[0];
rz(-0.11798141) q[0];
sx q[0];
rz(0.82103658) q[0];
rz(3.0631284) q[1];
sx q[1];
rz(-1.4469701) q[1];
sx q[1];
rz(-2.1441114) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8686912) q[0];
sx q[0];
rz(-1.2558381) q[0];
sx q[0];
rz(0.64906831) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1053379) q[2];
sx q[2];
rz(-1.6299575) q[2];
sx q[2];
rz(-1.349337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4892024) q[1];
sx q[1];
rz(-1.308131) q[1];
sx q[1];
rz(1.0242382) q[1];
rz(3.0313052) q[3];
sx q[3];
rz(-0.48088851) q[3];
sx q[3];
rz(0.65910027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0492101) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(-1.1537665) q[2];
rz(-2.9554101) q[3];
sx q[3];
rz(-0.081705339) q[3];
sx q[3];
rz(2.8103099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1368197) q[0];
sx q[0];
rz(-2.3240219) q[0];
sx q[0];
rz(1.9983043) q[0];
rz(-1.1413057) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(0.51796651) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9310936) q[0];
sx q[0];
rz(-0.8862952) q[0];
sx q[0];
rz(2.4263591) q[0];
rz(-pi) q[1];
rz(-0.00092351726) q[2];
sx q[2];
rz(-1.5557655) q[2];
sx q[2];
rz(1.3236486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8570366) q[1];
sx q[1];
rz(-1.9003881) q[1];
sx q[1];
rz(2.8524141) q[1];
rz(-2.6916396) q[3];
sx q[3];
rz(-2.8128689) q[3];
sx q[3];
rz(-2.4462819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(-2.8138568) q[2];
rz(1.1945126) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(-0.85668844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1389393) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(-2.5797381) q[0];
rz(1.6506763) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(3.0432826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.597991) q[0];
sx q[0];
rz(-2.7181135) q[0];
sx q[0];
rz(1.4844358) q[0];
rz(0.00064050015) q[2];
sx q[2];
rz(-1.570669) q[2];
sx q[2];
rz(-1.887111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6209539) q[1];
sx q[1];
rz(-2.6700182) q[1];
sx q[1];
rz(-3.0388799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0344072) q[3];
sx q[3];
rz(-1.1157728) q[3];
sx q[3];
rz(-0.95038271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0564698) q[2];
sx q[2];
rz(-2.9462908) q[2];
sx q[2];
rz(2.0758212) q[2];
rz(-0.39984518) q[3];
sx q[3];
rz(-2.6078434) q[3];
sx q[3];
rz(-1.2679509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9999076) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(-1.4578693) q[0];
rz(-2.0211925) q[1];
sx q[1];
rz(-0.1381865) q[1];
sx q[1];
rz(2.8021326) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65050478) q[0];
sx q[0];
rz(-1.5575214) q[0];
sx q[0];
rz(0.076939452) q[0];
rz(0.55468126) q[2];
sx q[2];
rz(-1.5590057) q[2];
sx q[2];
rz(-1.5654711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5875747) q[1];
sx q[1];
rz(-1.6563935) q[1];
sx q[1];
rz(3.1399591) q[1];
rz(0.27642997) q[3];
sx q[3];
rz(-0.35110858) q[3];
sx q[3];
rz(-1.029226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3642984) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(2.778229) q[2];
rz(-1.0904788) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(2.0232078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61984396) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(-2.2186665) q[0];
rz(1.482831) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(-3.1373851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7711346) q[0];
sx q[0];
rz(-0.97386375) q[0];
sx q[0];
rz(0.93907066) q[0];
rz(-pi) q[1];
rz(-1.9816859) q[2];
sx q[2];
rz(-1.5653603) q[2];
sx q[2];
rz(1.5559352) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9862822) q[1];
sx q[1];
rz(-2.661507) q[1];
sx q[1];
rz(-1.4309149) q[1];
rz(-1.5245304) q[3];
sx q[3];
rz(-2.1827201) q[3];
sx q[3];
rz(0.16976742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5620455) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(1.2008249) q[2];
rz(1.3935401) q[3];
sx q[3];
rz(-0.0037007185) q[3];
sx q[3];
rz(2.4414731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30018184) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(-1.2196983) q[0];
rz(1.8118743) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(-2.9776998) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3541479) q[0];
sx q[0];
rz(-1.3992926) q[0];
sx q[0];
rz(1.9308912) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8437949) q[2];
sx q[2];
rz(-2.1381209) q[2];
sx q[2];
rz(1.8253872) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6170609) q[1];
sx q[1];
rz(-1.5111898) q[1];
sx q[1];
rz(1.4527133) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71299841) q[3];
sx q[3];
rz(-2.2621691) q[3];
sx q[3];
rz(2.612243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4280052) q[2];
sx q[2];
rz(-3.0516477) q[2];
sx q[2];
rz(0.62291992) q[2];
rz(3.0981787) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(-2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49122214) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(-2.6553335) q[0];
rz(-2.4468415) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(0.4756701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9566475) q[0];
sx q[0];
rz(-1.534919) q[0];
sx q[0];
rz(0.033680276) q[0];
rz(-pi) q[1];
rz(-0.86464244) q[2];
sx q[2];
rz(-1.5448031) q[2];
sx q[2];
rz(-3.1033422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.039363843) q[1];
sx q[1];
rz(-2.1783531) q[1];
sx q[1];
rz(-1.5786912) q[1];
rz(-pi) q[2];
rz(2.6127897) q[3];
sx q[3];
rz(-2.7861054) q[3];
sx q[3];
rz(0.58099174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4645369) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(1.8439058) q[2];
rz(-1.8929947) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(-2.7738074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1260592) q[0];
sx q[0];
rz(-1.561469) q[0];
sx q[0];
rz(-1.4373686) q[0];
rz(0.88232782) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(-1.5346943) q[2];
sx q[2];
rz(-1.4786199) q[2];
sx q[2];
rz(0.27603966) q[2];
rz(-2.9024057) q[3];
sx q[3];
rz(-1.5755972) q[3];
sx q[3];
rz(1.5920873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
