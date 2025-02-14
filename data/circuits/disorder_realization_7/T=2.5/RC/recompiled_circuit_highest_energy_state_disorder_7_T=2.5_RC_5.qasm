OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.67686358) q[0];
sx q[0];
rz(-1.6114177) q[0];
sx q[0];
rz(-2.8443008) q[0];
rz(-2.3625506) q[1];
sx q[1];
rz(-2.9809597) q[1];
sx q[1];
rz(1.7056486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4592847) q[0];
sx q[0];
rz(-0.78792446) q[0];
sx q[0];
rz(0.81972229) q[0];
x q[1];
rz(1.3117242) q[2];
sx q[2];
rz(-1.6951778) q[2];
sx q[2];
rz(-0.69394116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85992009) q[1];
sx q[1];
rz(-1.8919196) q[1];
sx q[1];
rz(-2.4008277) q[1];
rz(-pi) q[2];
rz(0.30835704) q[3];
sx q[3];
rz(-1.3801818) q[3];
sx q[3];
rz(0.38577831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5433189) q[2];
sx q[2];
rz(-0.42793772) q[2];
sx q[2];
rz(-2.9738026) q[2];
rz(-1.1281891) q[3];
sx q[3];
rz(-1.5725243) q[3];
sx q[3];
rz(-1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0447277) q[0];
sx q[0];
rz(-2.4714578) q[0];
sx q[0];
rz(-2.0665533) q[0];
rz(1.3661522) q[1];
sx q[1];
rz(-1.7551883) q[1];
sx q[1];
rz(1.3137821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15720651) q[0];
sx q[0];
rz(-1.7922548) q[0];
sx q[0];
rz(-2.8048672) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4013675) q[2];
sx q[2];
rz(-1.3870948) q[2];
sx q[2];
rz(1.4445164) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35006501) q[1];
sx q[1];
rz(-2.413837) q[1];
sx q[1];
rz(1.6863053) q[1];
x q[2];
rz(1.1978549) q[3];
sx q[3];
rz(-0.35344703) q[3];
sx q[3];
rz(1.272205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65110597) q[2];
sx q[2];
rz(-1.3765114) q[2];
sx q[2];
rz(-2.7172078) q[2];
rz(-0.73244798) q[3];
sx q[3];
rz(-1.4832352) q[3];
sx q[3];
rz(1.6781767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3378147) q[0];
sx q[0];
rz(-2.4265899) q[0];
sx q[0];
rz(-2.6388229) q[0];
rz(0.94364199) q[1];
sx q[1];
rz(-2.4991401) q[1];
sx q[1];
rz(0.81519333) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749933) q[0];
sx q[0];
rz(-0.19375676) q[0];
sx q[0];
rz(2.8961477) q[0];
rz(-0.72702144) q[2];
sx q[2];
rz(-0.15910782) q[2];
sx q[2];
rz(-0.66691676) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.013061831) q[1];
sx q[1];
rz(-1.8327886) q[1];
sx q[1];
rz(0.2688516) q[1];
x q[2];
rz(0.56166537) q[3];
sx q[3];
rz(-0.9540073) q[3];
sx q[3];
rz(-1.3514618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0649071) q[2];
sx q[2];
rz(-1.6904597) q[2];
sx q[2];
rz(-1.1265075) q[2];
rz(1.3569776) q[3];
sx q[3];
rz(-0.49961909) q[3];
sx q[3];
rz(0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821238) q[0];
sx q[0];
rz(-0.72002763) q[0];
sx q[0];
rz(0.35991392) q[0];
rz(-2.9008046) q[1];
sx q[1];
rz(-1.9937932) q[1];
sx q[1];
rz(-0.37318939) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68615896) q[0];
sx q[0];
rz(-2.9669274) q[0];
sx q[0];
rz(-2.1771978) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56314205) q[2];
sx q[2];
rz(-1.3712607) q[2];
sx q[2];
rz(-0.89479252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9149941) q[1];
sx q[1];
rz(-2.0862596) q[1];
sx q[1];
rz(1.1771918) q[1];
rz(-pi) q[2];
rz(-0.4469965) q[3];
sx q[3];
rz(-1.3986386) q[3];
sx q[3];
rz(-0.76154532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.084426247) q[2];
sx q[2];
rz(-1.7095704) q[2];
sx q[2];
rz(0.90844321) q[2];
rz(-1.069979) q[3];
sx q[3];
rz(-1.1845651) q[3];
sx q[3];
rz(2.158304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.1835566) q[0];
sx q[0];
rz(-2.2479842) q[0];
sx q[0];
rz(0.91304427) q[0];
rz(-0.68963447) q[1];
sx q[1];
rz(-2.2328186) q[1];
sx q[1];
rz(2.3675809) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1415187) q[0];
sx q[0];
rz(-2.1545) q[0];
sx q[0];
rz(1.0144177) q[0];
rz(-pi) q[1];
x q[1];
rz(0.795587) q[2];
sx q[2];
rz(-1.5930111) q[2];
sx q[2];
rz(1.8088248) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23394984) q[1];
sx q[1];
rz(-1.7617258) q[1];
sx q[1];
rz(1.1348866) q[1];
x q[2];
rz(0.048581913) q[3];
sx q[3];
rz(-2.5932199) q[3];
sx q[3];
rz(-2.2480132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.245605) q[2];
sx q[2];
rz(-1.7481123) q[2];
sx q[2];
rz(-0.73949933) q[2];
rz(1.6351522) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(0.80772775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5189811) q[0];
sx q[0];
rz(-2.3464572) q[0];
sx q[0];
rz(-1.125289) q[0];
rz(3.0033424) q[1];
sx q[1];
rz(-1.467265) q[1];
sx q[1];
rz(1.5743871) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.137321) q[0];
sx q[0];
rz(-1.9372276) q[0];
sx q[0];
rz(-1.9204813) q[0];
rz(-pi) q[1];
rz(1.862941) q[2];
sx q[2];
rz(-1.9685479) q[2];
sx q[2];
rz(1.6145371) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9003657) q[1];
sx q[1];
rz(-1.8292184) q[1];
sx q[1];
rz(1.4421806) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7683783) q[3];
sx q[3];
rz(-2.1113331) q[3];
sx q[3];
rz(0.61321248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77341998) q[2];
sx q[2];
rz(-1.0680826) q[2];
sx q[2];
rz(-1.0991905) q[2];
rz(-2.1558971) q[3];
sx q[3];
rz(-1.516195) q[3];
sx q[3];
rz(-0.43258468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8428335) q[0];
sx q[0];
rz(-2.1465813) q[0];
sx q[0];
rz(0.80879912) q[0];
rz(0.66566268) q[1];
sx q[1];
rz(-0.97949615) q[1];
sx q[1];
rz(-2.1601802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8869276) q[0];
sx q[0];
rz(-1.6279632) q[0];
sx q[0];
rz(0.99941038) q[0];
rz(-pi) q[1];
rz(-1.8265884) q[2];
sx q[2];
rz(-1.0793387) q[2];
sx q[2];
rz(-2.0602399) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3700404) q[1];
sx q[1];
rz(-2.0345033) q[1];
sx q[1];
rz(0.15180219) q[1];
x q[2];
rz(0.1372758) q[3];
sx q[3];
rz(-0.90291728) q[3];
sx q[3];
rz(-1.5811416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.064528) q[2];
sx q[2];
rz(-1.6669824) q[2];
sx q[2];
rz(0.36901739) q[2];
rz(-1.4024233) q[3];
sx q[3];
rz(-0.88715059) q[3];
sx q[3];
rz(-1.3212475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8333261) q[0];
sx q[0];
rz(-2.3127191) q[0];
sx q[0];
rz(2.8114124) q[0];
rz(-0.45686832) q[1];
sx q[1];
rz(-2.7062682) q[1];
sx q[1];
rz(0.59828573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56494026) q[0];
sx q[0];
rz(-2.1974753) q[0];
sx q[0];
rz(-2.1184854) q[0];
rz(2.0399698) q[2];
sx q[2];
rz(-1.507505) q[2];
sx q[2];
rz(0.67055145) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48409325) q[1];
sx q[1];
rz(-1.7998905) q[1];
sx q[1];
rz(0.34219663) q[1];
x q[2];
rz(2.1586447) q[3];
sx q[3];
rz(-2.1390669) q[3];
sx q[3];
rz(-2.1592888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8699708) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(3.0769707) q[2];
rz(-1.6914852) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(1.5959285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86218631) q[0];
sx q[0];
rz(-2.8074844) q[0];
sx q[0];
rz(-2.620328) q[0];
rz(2.0817256) q[1];
sx q[1];
rz(-1.9797378) q[1];
sx q[1];
rz(-1.0313787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17612621) q[0];
sx q[0];
rz(-2.3076008) q[0];
sx q[0];
rz(1.6710207) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8487186) q[2];
sx q[2];
rz(-1.4559064) q[2];
sx q[2];
rz(0.024115901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4881812) q[1];
sx q[1];
rz(-1.0489042) q[1];
sx q[1];
rz(-0.33409358) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4724588) q[3];
sx q[3];
rz(-1.0957484) q[3];
sx q[3];
rz(-0.46365427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3970268) q[2];
sx q[2];
rz(-0.72781813) q[2];
sx q[2];
rz(-3.0542206) q[2];
rz(-1.2157863) q[3];
sx q[3];
rz(-1.4601424) q[3];
sx q[3];
rz(-0.17010918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4755197) q[0];
sx q[0];
rz(-0.89101321) q[0];
sx q[0];
rz(-2.0024306) q[0];
rz(-0.49925223) q[1];
sx q[1];
rz(-1.6923994) q[1];
sx q[1];
rz(1.8081236) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069786929) q[0];
sx q[0];
rz(-0.13955586) q[0];
sx q[0];
rz(-2.4371171) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3158445) q[2];
sx q[2];
rz(-1.5226018) q[2];
sx q[2];
rz(2.9914837) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6302609) q[1];
sx q[1];
rz(-1.7862661) q[1];
sx q[1];
rz(-2.8829619) q[1];
rz(-pi) q[2];
rz(0.74204294) q[3];
sx q[3];
rz(-2.2931406) q[3];
sx q[3];
rz(2.6076041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6405876) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(2.8583543) q[2];
rz(1.2011652) q[3];
sx q[3];
rz(-2.4656506) q[3];
sx q[3];
rz(1.0454319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129811) q[0];
sx q[0];
rz(-2.3647478) q[0];
sx q[0];
rz(1.1051529) q[0];
rz(-0.10440566) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(-1.3731643) q[2];
sx q[2];
rz(-2.0680972) q[2];
sx q[2];
rz(1.8152449) q[2];
rz(2.4734409) q[3];
sx q[3];
rz(-1.9037608) q[3];
sx q[3];
rz(-1.9477061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
