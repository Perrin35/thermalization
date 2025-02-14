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
rz(-1.4359441) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4592847) q[0];
sx q[0];
rz(-0.78792446) q[0];
sx q[0];
rz(-2.3218704) q[0];
x q[1];
rz(-3.0129635) q[2];
sx q[2];
rz(-1.8278215) q[2];
sx q[2];
rz(0.84398735) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7117225) q[1];
sx q[1];
rz(-0.87583576) q[1];
sx q[1];
rz(1.9942787) q[1];
rz(-pi) q[2];
rz(-2.5752742) q[3];
sx q[3];
rz(-2.7806823) q[3];
sx q[3];
rz(-2.4931815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5433189) q[2];
sx q[2];
rz(-2.7136549) q[2];
sx q[2];
rz(-2.9738026) q[2];
rz(-2.0134036) q[3];
sx q[3];
rz(-1.5725243) q[3];
sx q[3];
rz(1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0447277) q[0];
sx q[0];
rz(-2.4714578) q[0];
sx q[0];
rz(-2.0665533) q[0];
rz(-1.7754405) q[1];
sx q[1];
rz(-1.7551883) q[1];
sx q[1];
rz(1.3137821) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15720651) q[0];
sx q[0];
rz(-1.7922548) q[0];
sx q[0];
rz(0.33672543) q[0];
rz(2.4046217) q[2];
sx q[2];
rz(-0.24925512) q[2];
sx q[2];
rz(0.94446206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35006501) q[1];
sx q[1];
rz(-0.72775562) q[1];
sx q[1];
rz(1.4552874) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9437378) q[3];
sx q[3];
rz(-0.35344703) q[3];
sx q[3];
rz(1.8693876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65110597) q[2];
sx q[2];
rz(-1.7650812) q[2];
sx q[2];
rz(-0.42438486) q[2];
rz(0.73244798) q[3];
sx q[3];
rz(-1.4832352) q[3];
sx q[3];
rz(-1.6781767) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80377793) q[0];
sx q[0];
rz(-0.71500272) q[0];
sx q[0];
rz(2.6388229) q[0];
rz(2.1979507) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(-2.3263993) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8250587) q[0];
sx q[0];
rz(-1.382917) q[0];
sx q[0];
rz(-1.618439) q[0];
x q[1];
rz(-1.6770467) q[2];
sx q[2];
rz(-1.4521404) q[2];
sx q[2];
rz(-0.066421631) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1285308) q[1];
sx q[1];
rz(-1.3088041) q[1];
sx q[1];
rz(2.872741) q[1];
x q[2];
rz(-2.2150119) q[3];
sx q[3];
rz(-0.80873064) q[3];
sx q[3];
rz(0.96264983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.076685585) q[2];
sx q[2];
rz(-1.6904597) q[2];
sx q[2];
rz(1.1265075) q[2];
rz(-1.7846151) q[3];
sx q[3];
rz(-2.6419736) q[3];
sx q[3];
rz(-0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821238) q[0];
sx q[0];
rz(-0.72002763) q[0];
sx q[0];
rz(-2.7816787) q[0];
rz(2.9008046) q[1];
sx q[1];
rz(-1.1477995) q[1];
sx q[1];
rz(-0.37318939) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8418331) q[0];
sx q[0];
rz(-1.4275121) q[0];
sx q[0];
rz(-0.10023198) q[0];
x q[1];
rz(2.5784506) q[2];
sx q[2];
rz(-1.7703319) q[2];
sx q[2];
rz(2.2468001) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9149941) q[1];
sx q[1];
rz(-1.0553331) q[1];
sx q[1];
rz(-1.1771918) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7612824) q[3];
sx q[3];
rz(-2.0107186) q[3];
sx q[3];
rz(-0.89118496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0571664) q[2];
sx q[2];
rz(-1.7095704) q[2];
sx q[2];
rz(0.90844321) q[2];
rz(1.069979) q[3];
sx q[3];
rz(-1.1845651) q[3];
sx q[3];
rz(-2.158304) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95803607) q[0];
sx q[0];
rz(-0.89360845) q[0];
sx q[0];
rz(2.2285484) q[0];
rz(2.4519582) q[1];
sx q[1];
rz(-2.2328186) q[1];
sx q[1];
rz(-0.77401179) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90092313) q[0];
sx q[0];
rz(-1.114448) q[0];
sx q[0];
rz(-2.4805446) q[0];
x q[1];
rz(0.031096259) q[2];
sx q[2];
rz(-0.79582873) q[2];
sx q[2];
rz(2.9253256) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8929114) q[1];
sx q[1];
rz(-1.9982583) q[1];
sx q[1];
rz(-0.21007725) q[1];
x q[2];
rz(-1.6004531) q[3];
sx q[3];
rz(-2.1184485) q[3];
sx q[3];
rz(0.95049196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.245605) q[2];
sx q[2];
rz(-1.3934803) q[2];
sx q[2];
rz(2.4020933) q[2];
rz(1.6351522) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(0.80772775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6226115) q[0];
sx q[0];
rz(-0.7951355) q[0];
sx q[0];
rz(-1.125289) q[0];
rz(-3.0033424) q[1];
sx q[1];
rz(-1.467265) q[1];
sx q[1];
rz(-1.5743871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7985021) q[0];
sx q[0];
rz(-2.6406086) q[0];
sx q[0];
rz(-0.72879099) q[0];
rz(-0.41344686) q[2];
sx q[2];
rz(-1.3020497) q[2];
sx q[2];
rz(-0.15971249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77273276) q[1];
sx q[1];
rz(-0.28801685) q[1];
sx q[1];
rz(-0.45175987) q[1];
rz(-pi) q[2];
rz(0.31611021) q[3];
sx q[3];
rz(-0.5721285) q[3];
sx q[3];
rz(2.8993949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3681727) q[2];
sx q[2];
rz(-1.0680826) q[2];
sx q[2];
rz(1.0991905) q[2];
rz(2.1558971) q[3];
sx q[3];
rz(-1.6253977) q[3];
sx q[3];
rz(-0.43258468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.2987591) q[0];
sx q[0];
rz(-2.1465813) q[0];
sx q[0];
rz(-0.80879912) q[0];
rz(0.66566268) q[1];
sx q[1];
rz(-0.97949615) q[1];
sx q[1];
rz(-2.1601802) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4047667) q[0];
sx q[0];
rz(-0.5739218) q[0];
sx q[0];
rz(-1.6762275) q[0];
x q[1];
rz(-2.6362474) q[2];
sx q[2];
rz(-1.7957558) q[2];
sx q[2];
rz(-2.7749429) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8725207) q[1];
sx q[1];
rz(-1.4351294) q[1];
sx q[1];
rz(2.0391463) q[1];
rz(-pi) q[2];
rz(3.0043169) q[3];
sx q[3];
rz(-0.90291728) q[3];
sx q[3];
rz(1.5811416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0770646) q[2];
sx q[2];
rz(-1.6669824) q[2];
sx q[2];
rz(-0.36901739) q[2];
rz(-1.7391694) q[3];
sx q[3];
rz(-2.2544421) q[3];
sx q[3];
rz(-1.3212475) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8333261) q[0];
sx q[0];
rz(-2.3127191) q[0];
sx q[0];
rz(-2.8114124) q[0];
rz(2.6847243) q[1];
sx q[1];
rz(-0.43532443) q[1];
sx q[1];
rz(2.5433069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56494026) q[0];
sx q[0];
rz(-0.9441174) q[0];
sx q[0];
rz(2.1184854) q[0];
rz(-pi) q[1];
rz(0.070934709) q[2];
sx q[2];
rz(-2.0389551) q[2];
sx q[2];
rz(2.2092961) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51914267) q[1];
sx q[1];
rz(-0.40928091) q[1];
sx q[1];
rz(0.60731387) q[1];
rz(-pi) q[2];
rz(-2.4264494) q[3];
sx q[3];
rz(-0.79350457) q[3];
sx q[3];
rz(1.8737829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2716219) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(-3.0769707) q[2];
rz(1.4501075) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(-1.5456642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86218631) q[0];
sx q[0];
rz(-2.8074844) q[0];
sx q[0];
rz(-0.52126467) q[0];
rz(2.0817256) q[1];
sx q[1];
rz(-1.9797378) q[1];
sx q[1];
rz(2.1102139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3272031) q[0];
sx q[0];
rz(-1.6449682) q[0];
sx q[0];
rz(0.73930862) q[0];
x q[1];
rz(2.7613369) q[2];
sx q[2];
rz(-0.31399841) q[2];
sx q[2];
rz(1.2316201) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0887275) q[1];
sx q[1];
rz(-1.8590312) q[1];
sx q[1];
rz(1.0239787) q[1];
rz(-pi) q[2];
rz(-2.9529753) q[3];
sx q[3];
rz(-0.48435703) q[3];
sx q[3];
rz(2.8903807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7445659) q[2];
sx q[2];
rz(-2.4137745) q[2];
sx q[2];
rz(-3.0542206) q[2];
rz(-1.2157863) q[3];
sx q[3];
rz(-1.6814503) q[3];
sx q[3];
rz(-2.9714835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66607296) q[0];
sx q[0];
rz(-2.2505794) q[0];
sx q[0];
rz(-1.1391621) q[0];
rz(-0.49925223) q[1];
sx q[1];
rz(-1.4491932) q[1];
sx q[1];
rz(1.3334691) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0718057) q[0];
sx q[0];
rz(-0.13955586) q[0];
sx q[0];
rz(-2.4371171) q[0];
rz(-pi) q[1];
rz(-1.3818327) q[2];
sx q[2];
rz(-0.25936959) q[2];
sx q[2];
rz(1.9036906) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51133174) q[1];
sx q[1];
rz(-1.7862661) q[1];
sx q[1];
rz(2.8829619) q[1];
rz(2.3995497) q[3];
sx q[3];
rz(-0.84845209) q[3];
sx q[3];
rz(-0.53398856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.501005) q[2];
sx q[2];
rz(-1.86684) q[2];
sx q[2];
rz(2.8583543) q[2];
rz(-1.9404274) q[3];
sx q[3];
rz(-2.4656506) q[3];
sx q[3];
rz(-2.0961608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.82861154) q[0];
sx q[0];
rz(-2.3647478) q[0];
sx q[0];
rz(1.1051529) q[0];
rz(-0.10440566) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(0.5055867) q[2];
sx q[2];
rz(-1.3973631) q[2];
sx q[2];
rz(0.33968795) q[2];
rz(1.9857882) q[3];
sx q[3];
rz(-2.1963051) q[3];
sx q[3];
rz(2.5122234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
