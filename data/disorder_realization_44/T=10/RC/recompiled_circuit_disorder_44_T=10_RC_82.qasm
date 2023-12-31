OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(3.0043998) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44088988) q[0];
sx q[0];
rz(-0.39632495) q[0];
sx q[0];
rz(-1.2896145) q[0];
rz(-pi) q[1];
rz(2.2719703) q[2];
sx q[2];
rz(-2.6343971) q[2];
sx q[2];
rz(1.6091572) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4094028) q[1];
sx q[1];
rz(-0.70530546) q[1];
sx q[1];
rz(2.3975055) q[1];
rz(-pi) q[2];
rz(2.8785273) q[3];
sx q[3];
rz(-2.3657626) q[3];
sx q[3];
rz(-2.8924931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(1.6254788) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(1.6158993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34823725) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(-3.120378) q[0];
rz(-1.1938098) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-0.83591998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0854557) q[0];
sx q[0];
rz(-1.5888927) q[0];
sx q[0];
rz(-1.7096814) q[0];
rz(2.1875728) q[2];
sx q[2];
rz(-0.6435794) q[2];
sx q[2];
rz(-1.8804903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0554725) q[1];
sx q[1];
rz(-1.8430084) q[1];
sx q[1];
rz(-2.230456) q[1];
rz(1.9637945) q[3];
sx q[3];
rz(-2.2713695) q[3];
sx q[3];
rz(-0.29004471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2772284) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(1.7956087) q[2];
rz(-2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(1.0536449) q[0];
rz(-1.9127649) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(0.4371117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8374098) q[0];
sx q[0];
rz(-1.279631) q[0];
sx q[0];
rz(-3.0512179) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0929298) q[2];
sx q[2];
rz(-0.46233593) q[2];
sx q[2];
rz(-2.1585652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8902953) q[1];
sx q[1];
rz(-1.1493059) q[1];
sx q[1];
rz(1.4189659) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4594853) q[3];
sx q[3];
rz(-2.6442332) q[3];
sx q[3];
rz(-1.5670083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.019471021) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(-1.0220698) q[2];
rz(1.2381037) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(-2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(0.98130256) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-0.19128004) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6879038) q[0];
sx q[0];
rz(-0.17246248) q[0];
sx q[0];
rz(1.0426636) q[0];
rz(-0.83696604) q[2];
sx q[2];
rz(-1.776473) q[2];
sx q[2];
rz(-0.93081805) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3595708) q[1];
sx q[1];
rz(-1.6732209) q[1];
sx q[1];
rz(2.3624079) q[1];
rz(0.55320923) q[3];
sx q[3];
rz(-2.534453) q[3];
sx q[3];
rz(0.34724423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(1.0106687) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8028832) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(0.55661911) q[0];
rz(-0.11511766) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-2.1690878) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33362493) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(-1.6385965) q[0];
rz(1.2178671) q[2];
sx q[2];
rz(-2.3595516) q[2];
sx q[2];
rz(-1.023934) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79341187) q[1];
sx q[1];
rz(-0.23985292) q[1];
sx q[1];
rz(-0.94437771) q[1];
x q[2];
rz(-1.6226107) q[3];
sx q[3];
rz(-2.5250146) q[3];
sx q[3];
rz(0.50293621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82289034) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.5931607) q[2];
rz(-1.7758153) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(-0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.71972972) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(-1.4259889) q[0];
rz(-1.0643719) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(-2.7672966) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4842589) q[0];
sx q[0];
rz(-1.9801635) q[0];
sx q[0];
rz(1.9166458) q[0];
x q[1];
rz(-0.46005581) q[2];
sx q[2];
rz(-0.54121491) q[2];
sx q[2];
rz(-0.40749007) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4243851) q[1];
sx q[1];
rz(-1.8339515) q[1];
sx q[1];
rz(2.5724263) q[1];
x q[2];
rz(-2.2603287) q[3];
sx q[3];
rz(-2.1173819) q[3];
sx q[3];
rz(1.5342086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8950243) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(-0.22917497) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(-2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36528698) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(-2.2348485) q[0];
rz(-2.0523741) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(1.3100756) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10446564) q[0];
sx q[0];
rz(-1.5195527) q[0];
sx q[0];
rz(-0.38802223) q[0];
rz(-1.2398948) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(2.9943525) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20565614) q[1];
sx q[1];
rz(-1.9034791) q[1];
sx q[1];
rz(-1.3693621) q[1];
x q[2];
rz(-0.044192627) q[3];
sx q[3];
rz(-0.5558388) q[3];
sx q[3];
rz(2.1293872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2157796) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(-2.8619134) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(-0.057597615) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16335547) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(-2.9220007) q[0];
rz(2.638468) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(0.84987744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06892878) q[0];
sx q[0];
rz(-1.7964296) q[0];
sx q[0];
rz(0.16107852) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4555898) q[2];
sx q[2];
rz(-0.68636471) q[2];
sx q[2];
rz(2.3442868) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0440158) q[1];
sx q[1];
rz(-2.3313287) q[1];
sx q[1];
rz(-0.047659831) q[1];
x q[2];
rz(1.1342808) q[3];
sx q[3];
rz(-0.28537649) q[3];
sx q[3];
rz(0.34261045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59051096) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.3809416) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6324156) q[0];
sx q[0];
rz(-0.89288765) q[0];
sx q[0];
rz(0.40503043) q[0];
rz(0.45267725) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(-1.8639494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687846) q[0];
sx q[0];
rz(-2.0622258) q[0];
sx q[0];
rz(1.8322893) q[0];
rz(2.4503166) q[2];
sx q[2];
rz(-1.9553767) q[2];
sx q[2];
rz(-2.861475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7840522) q[1];
sx q[1];
rz(-1.6077542) q[1];
sx q[1];
rz(2.0164911) q[1];
rz(-pi) q[2];
rz(-2.7550335) q[3];
sx q[3];
rz(-2.2694526) q[3];
sx q[3];
rz(-0.78702918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8032802) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(2.7837616) q[2];
rz(1.4194277) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2232067) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(0.11225587) q[0];
rz(-0.90011251) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(-0.21044883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4608599) q[0];
sx q[0];
rz(-2.0740777) q[0];
sx q[0];
rz(-1.7778648) q[0];
rz(-pi) q[1];
rz(-1.9824355) q[2];
sx q[2];
rz(-1.4201846) q[2];
sx q[2];
rz(2.6580236) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3810972) q[1];
sx q[1];
rz(-1.7777182) q[1];
sx q[1];
rz(2.8269672) q[1];
x q[2];
rz(-2.8912192) q[3];
sx q[3];
rz(-1.8009406) q[3];
sx q[3];
rz(1.9104513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9615053) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(1.2735584) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(-2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56959854) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(2.3251484) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(1.1007166) q[2];
sx q[2];
rz(-1.7190949) q[2];
sx q[2];
rz(0.20863056) q[2];
rz(2.7995085) q[3];
sx q[3];
rz(-2.2707006) q[3];
sx q[3];
rz(-2.0778098) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
