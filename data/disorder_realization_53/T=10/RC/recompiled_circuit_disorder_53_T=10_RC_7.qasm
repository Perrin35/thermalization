OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18908137) q[0];
sx q[0];
rz(-0.12726769) q[0];
sx q[0];
rz(-2.1531818) q[0];
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(1.3487863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.792359) q[0];
sx q[0];
rz(-1.4548737) q[0];
sx q[0];
rz(0.16616343) q[0];
x q[1];
rz(1.4853391) q[2];
sx q[2];
rz(-1.8329617) q[2];
sx q[2];
rz(-2.0567577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3218282) q[1];
sx q[1];
rz(-1.0294224) q[1];
sx q[1];
rz(0.86408792) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80356055) q[3];
sx q[3];
rz(-0.27419146) q[3];
sx q[3];
rz(2.1823332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(-2.7963426) q[2];
rz(-2.9521862) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(-0.96864831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(0.35821113) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(-2.125724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43243877) q[0];
sx q[0];
rz(-1.0854939) q[0];
sx q[0];
rz(-2.0251459) q[0];
x q[1];
rz(-2.303896) q[2];
sx q[2];
rz(-2.838755) q[2];
sx q[2];
rz(0.31485117) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3381172) q[1];
sx q[1];
rz(-2.7403767) q[1];
sx q[1];
rz(-0.34715279) q[1];
rz(-pi) q[2];
rz(2.0738828) q[3];
sx q[3];
rz(-1.9670608) q[3];
sx q[3];
rz(2.577142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7887855) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-2.0324198) q[2];
rz(0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(-3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1576841) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(-2.2614959) q[0];
rz(1.8799211) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(2.9601011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49849579) q[0];
sx q[0];
rz(-2.4884014) q[0];
sx q[0];
rz(0.80723395) q[0];
rz(-pi) q[1];
rz(0.7736189) q[2];
sx q[2];
rz(-2.3256989) q[2];
sx q[2];
rz(1.1676163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.03141244) q[1];
sx q[1];
rz(-1.8437779) q[1];
sx q[1];
rz(-0.024755342) q[1];
x q[2];
rz(0.39748945) q[3];
sx q[3];
rz(-2.6263413) q[3];
sx q[3];
rz(1.6911563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.024753831) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(1.3428358) q[2];
rz(1.8163619) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(2.0480115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745859) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(0.68853199) q[0];
rz(-0.21903285) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(-2.9825488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6200136) q[0];
sx q[0];
rz(-2.0984681) q[0];
sx q[0];
rz(2.2844727) q[0];
rz(-pi) q[1];
rz(0.7775457) q[2];
sx q[2];
rz(-0.50707196) q[2];
sx q[2];
rz(2.5271202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98475458) q[1];
sx q[1];
rz(-1.5459237) q[1];
sx q[1];
rz(-1.7825885) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3499447) q[3];
sx q[3];
rz(-2.0851118) q[3];
sx q[3];
rz(0.056003464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8461385) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(0.42567483) q[2];
rz(-0.019429026) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(-0.69452906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6762125) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(-0.12839578) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(-2.593186) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5760561) q[0];
sx q[0];
rz(-2.5649374) q[0];
sx q[0];
rz(-1.499349) q[0];
rz(-pi) q[1];
rz(1.7565281) q[2];
sx q[2];
rz(-0.81132946) q[2];
sx q[2];
rz(3.082049) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4244528) q[1];
sx q[1];
rz(-1.6315178) q[1];
sx q[1];
rz(2.0408003) q[1];
rz(-3.0382285) q[3];
sx q[3];
rz(-1.169239) q[3];
sx q[3];
rz(-0.38927024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8473062) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(1.4667286) q[2];
rz(2.0560125) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(-0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(2.0423245) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(-0.99463314) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69840616) q[0];
sx q[0];
rz(-0.61265677) q[0];
sx q[0];
rz(-1.1820656) q[0];
x q[1];
rz(0.0066130916) q[2];
sx q[2];
rz(-0.50351876) q[2];
sx q[2];
rz(2.7343482) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99018807) q[1];
sx q[1];
rz(-1.2328887) q[1];
sx q[1];
rz(0.0068454725) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1676222) q[3];
sx q[3];
rz(-2.2304428) q[3];
sx q[3];
rz(-2.6350104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(-0.4449521) q[2];
rz(2.3343202) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(-0.81594938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31535661) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(-0.89609599) q[0];
rz(0.90944666) q[1];
sx q[1];
rz(-2.8912631) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5298115) q[0];
sx q[0];
rz(-1.6336332) q[0];
sx q[0];
rz(-1.9301374) q[0];
rz(-1.0117815) q[2];
sx q[2];
rz(-0.85306877) q[2];
sx q[2];
rz(0.43874028) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8304886) q[1];
sx q[1];
rz(-1.6357592) q[1];
sx q[1];
rz(0.91598367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92618561) q[3];
sx q[3];
rz(-1.965596) q[3];
sx q[3];
rz(1.3046164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(1.9141076) q[2];
rz(-3.098439) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(2.9149122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9880923) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(-0.41326997) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(0.73910284) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9773213) q[0];
sx q[0];
rz(-0.56557206) q[0];
sx q[0];
rz(-0.092202734) q[0];
rz(-1.5831645) q[2];
sx q[2];
rz(-2.6855199) q[2];
sx q[2];
rz(-0.2218483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6369789) q[1];
sx q[1];
rz(-2.4006872) q[1];
sx q[1];
rz(1.3330589) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8352175) q[3];
sx q[3];
rz(-1.6248871) q[3];
sx q[3];
rz(-0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3070613) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(2.0397662) q[2];
rz(0.39673355) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(-0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48801625) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(2.5226412) q[0];
rz(-2.1221819) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(-2.6143262) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70452481) q[0];
sx q[0];
rz(-1.8509812) q[0];
sx q[0];
rz(1.2824476) q[0];
x q[1];
rz(-0.68248827) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(1.3860821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.077721715) q[1];
sx q[1];
rz(-1.9943024) q[1];
sx q[1];
rz(-1.0788171) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91068565) q[3];
sx q[3];
rz(-1.3798957) q[3];
sx q[3];
rz(-1.1235352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8573389) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(-2.3305317) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60944027) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(2.6623181) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-2.9634109) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44706599) q[0];
sx q[0];
rz(-0.94213001) q[0];
sx q[0];
rz(0.64283128) q[0];
x q[1];
rz(-0.41583305) q[2];
sx q[2];
rz(-0.73199474) q[2];
sx q[2];
rz(0.40643613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61160117) q[1];
sx q[1];
rz(-1.0582663) q[1];
sx q[1];
rz(-0.26228735) q[1];
rz(1.6931157) q[3];
sx q[3];
rz(-1.1432075) q[3];
sx q[3];
rz(-1.3631671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46897727) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(-0.82328063) q[2];
rz(3.0205884) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2330033) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(1.9227149) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(-3.1115816) q[2];
sx q[2];
rz(-1.9297615) q[2];
sx q[2];
rz(2.161138) q[2];
rz(2.0291438) q[3];
sx q[3];
rz(-2.4670798) q[3];
sx q[3];
rz(-2.3453875) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
