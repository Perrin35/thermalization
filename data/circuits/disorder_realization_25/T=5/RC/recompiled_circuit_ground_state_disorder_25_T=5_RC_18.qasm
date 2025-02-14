OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8018262) q[0];
sx q[0];
rz(-2.5751994) q[0];
sx q[0];
rz(1.9422148) q[0];
rz(-0.83130032) q[1];
sx q[1];
rz(-1.5273233) q[1];
sx q[1];
rz(-0.49973127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2466149) q[0];
sx q[0];
rz(-2.3736989) q[0];
sx q[0];
rz(1.8706066) q[0];
rz(-2.9414177) q[2];
sx q[2];
rz(-1.9250637) q[2];
sx q[2];
rz(-1.406671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.079165212) q[1];
sx q[1];
rz(-1.2794331) q[1];
sx q[1];
rz(0.77701648) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2890782) q[3];
sx q[3];
rz(-1.6073462) q[3];
sx q[3];
rz(0.73383777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97303331) q[2];
sx q[2];
rz(-1.7567822) q[2];
sx q[2];
rz(0.92156571) q[2];
rz(-1.5714931) q[3];
sx q[3];
rz(-0.67155963) q[3];
sx q[3];
rz(-2.4422372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2498995) q[0];
sx q[0];
rz(-1.0279011) q[0];
sx q[0];
rz(1.4693042) q[0];
rz(1.4768614) q[1];
sx q[1];
rz(-1.7988484) q[1];
sx q[1];
rz(-2.6254168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75137532) q[0];
sx q[0];
rz(-1.1773407) q[0];
sx q[0];
rz(-0.1050001) q[0];
rz(-1.2391813) q[2];
sx q[2];
rz(-1.3388763) q[2];
sx q[2];
rz(0.9882016) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.33715926) q[1];
sx q[1];
rz(-0.26072956) q[1];
sx q[1];
rz(-2.8689137) q[1];
x q[2];
rz(-2.6593239) q[3];
sx q[3];
rz(-2.2029556) q[3];
sx q[3];
rz(-1.9763115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76253647) q[2];
sx q[2];
rz(-1.2006589) q[2];
sx q[2];
rz(2.8442247) q[2];
rz(0.51658806) q[3];
sx q[3];
rz(-1.1597495) q[3];
sx q[3];
rz(0.62492257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49279889) q[0];
sx q[0];
rz(-0.53557098) q[0];
sx q[0];
rz(0.18185644) q[0];
rz(-2.6975373) q[1];
sx q[1];
rz(-2.2477138) q[1];
sx q[1];
rz(3.0603538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59653283) q[0];
sx q[0];
rz(-1.1756011) q[0];
sx q[0];
rz(-1.5395274) q[0];
rz(-1.2368579) q[2];
sx q[2];
rz(-2.4231964) q[2];
sx q[2];
rz(-0.12018724) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4887152) q[1];
sx q[1];
rz(-1.3849568) q[1];
sx q[1];
rz(-1.0936827) q[1];
rz(1.3879561) q[3];
sx q[3];
rz(-0.5683848) q[3];
sx q[3];
rz(0.050384132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2150725) q[2];
sx q[2];
rz(-0.36538616) q[2];
sx q[2];
rz(-0.60058769) q[2];
rz(1.8961204) q[3];
sx q[3];
rz(-0.54491091) q[3];
sx q[3];
rz(-0.039552461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1323701) q[0];
sx q[0];
rz(-2.4339269) q[0];
sx q[0];
rz(0.02455499) q[0];
rz(-2.764616) q[1];
sx q[1];
rz(-1.8507277) q[1];
sx q[1];
rz(2.7797508) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0252286) q[0];
sx q[0];
rz(-1.4170237) q[0];
sx q[0];
rz(-0.31391252) q[0];
rz(-pi) q[1];
rz(-0.078388647) q[2];
sx q[2];
rz(-1.196821) q[2];
sx q[2];
rz(-1.0678991) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3256073) q[1];
sx q[1];
rz(-2.4802842) q[1];
sx q[1];
rz(-0.2925175) q[1];
x q[2];
rz(-0.41976069) q[3];
sx q[3];
rz(-1.1294522) q[3];
sx q[3];
rz(3.0361058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5234066) q[2];
sx q[2];
rz(-1.2647102) q[2];
sx q[2];
rz(-0.82320172) q[2];
rz(-0.19436714) q[3];
sx q[3];
rz(-0.40209642) q[3];
sx q[3];
rz(-0.91485867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.119656) q[0];
sx q[0];
rz(-1.5384262) q[0];
sx q[0];
rz(-0.62171474) q[0];
rz(-0.74553982) q[1];
sx q[1];
rz(-2.3108683) q[1];
sx q[1];
rz(-3.0527557) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1232077) q[0];
sx q[0];
rz(-1.443571) q[0];
sx q[0];
rz(1.6379665) q[0];
rz(-2.2584469) q[2];
sx q[2];
rz(-2.7252203) q[2];
sx q[2];
rz(2.1379545) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4712661) q[1];
sx q[1];
rz(-0.85976344) q[1];
sx q[1];
rz(1.9416351) q[1];
rz(-pi) q[2];
rz(-2.9778175) q[3];
sx q[3];
rz(-1.9036674) q[3];
sx q[3];
rz(0.44595697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2290153) q[2];
sx q[2];
rz(-1.6385767) q[2];
sx q[2];
rz(-1.5360606) q[2];
rz(2.33365) q[3];
sx q[3];
rz(-2.2980821) q[3];
sx q[3];
rz(1.1725175) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87869969) q[0];
sx q[0];
rz(-1.7182588) q[0];
sx q[0];
rz(0.95243564) q[0];
rz(1.7772504) q[1];
sx q[1];
rz(-1.2063426) q[1];
sx q[1];
rz(0.84699026) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5085707) q[0];
sx q[0];
rz(-1.428033) q[0];
sx q[0];
rz(1.6667913) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2405841) q[2];
sx q[2];
rz(-0.80494138) q[2];
sx q[2];
rz(3.0778468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0274732) q[1];
sx q[1];
rz(-0.67016685) q[1];
sx q[1];
rz(0.7188188) q[1];
rz(-1.5853508) q[3];
sx q[3];
rz(-1.4676009) q[3];
sx q[3];
rz(-1.1905014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9534037) q[2];
sx q[2];
rz(-1.0651411) q[2];
sx q[2];
rz(0.38189253) q[2];
rz(1.1000819) q[3];
sx q[3];
rz(-0.81172687) q[3];
sx q[3];
rz(-2.7303117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5683658) q[0];
sx q[0];
rz(-1.8170284) q[0];
sx q[0];
rz(2.6439917) q[0];
rz(0.35500232) q[1];
sx q[1];
rz(-1.4811131) q[1];
sx q[1];
rz(-2.3752046) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9123133) q[0];
sx q[0];
rz(-1.6054432) q[0];
sx q[0];
rz(0.51280068) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1670157) q[2];
sx q[2];
rz(-1.2082074) q[2];
sx q[2];
rz(0.18747839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38098225) q[1];
sx q[1];
rz(-2.8462324) q[1];
sx q[1];
rz(2.9350314) q[1];
x q[2];
rz(2.4637632) q[3];
sx q[3];
rz(-1.3732217) q[3];
sx q[3];
rz(3.1238926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7327205) q[2];
sx q[2];
rz(-1.5343821) q[2];
sx q[2];
rz(1.4917779) q[2];
rz(-2.0626227) q[3];
sx q[3];
rz(-1.8337199) q[3];
sx q[3];
rz(2.5274247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.806818) q[0];
sx q[0];
rz(-0.65009999) q[0];
sx q[0];
rz(-2.7244869) q[0];
rz(-0.053622309) q[1];
sx q[1];
rz(-1.5258748) q[1];
sx q[1];
rz(2.0448304) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7863207) q[0];
sx q[0];
rz(-1.6910166) q[0];
sx q[0];
rz(-0.97697301) q[0];
rz(-pi) q[1];
rz(-1.4545733) q[2];
sx q[2];
rz(-2.6781278) q[2];
sx q[2];
rz(1.506724) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.027710304) q[1];
sx q[1];
rz(-1.405763) q[1];
sx q[1];
rz(1.5636958) q[1];
x q[2];
rz(-0.18713672) q[3];
sx q[3];
rz(-0.68681648) q[3];
sx q[3];
rz(0.88469626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2199478) q[2];
sx q[2];
rz(-1.767445) q[2];
sx q[2];
rz(-2.5448223) q[2];
rz(2.2611332) q[3];
sx q[3];
rz(-1.7170649) q[3];
sx q[3];
rz(-0.60976353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6841458) q[0];
sx q[0];
rz(-0.54231751) q[0];
sx q[0];
rz(2.82161) q[0];
rz(2.7931702) q[1];
sx q[1];
rz(-2.6967144) q[1];
sx q[1];
rz(-0.50331032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5849294) q[0];
sx q[0];
rz(-1.4787424) q[0];
sx q[0];
rz(-1.7199055) q[0];
rz(-1.032802) q[2];
sx q[2];
rz(-0.82394407) q[2];
sx q[2];
rz(1.2822373) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.383068) q[1];
sx q[1];
rz(-0.43979859) q[1];
sx q[1];
rz(1.7563933) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4473404) q[3];
sx q[3];
rz(-2.0887825) q[3];
sx q[3];
rz(2.4855465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2549501) q[2];
sx q[2];
rz(-0.4526259) q[2];
sx q[2];
rz(-2.5060999) q[2];
rz(-1.6057711) q[3];
sx q[3];
rz(-1.7255892) q[3];
sx q[3];
rz(-3.0543069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4601634) q[0];
sx q[0];
rz(-2.5458113) q[0];
sx q[0];
rz(-1.7728565) q[0];
rz(0.5303371) q[1];
sx q[1];
rz(-1.4884721) q[1];
sx q[1];
rz(2.368685) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.146713) q[0];
sx q[0];
rz(-1.9582703) q[0];
sx q[0];
rz(-2.4058624) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10112986) q[2];
sx q[2];
rz(-0.79532184) q[2];
sx q[2];
rz(2.1908608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1369644) q[1];
sx q[1];
rz(-1.1234979) q[1];
sx q[1];
rz(-2.9699259) q[1];
x q[2];
rz(1.1201376) q[3];
sx q[3];
rz(-1.3067596) q[3];
sx q[3];
rz(-0.23917374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7698001) q[2];
sx q[2];
rz(-0.95855203) q[2];
sx q[2];
rz(-0.27251631) q[2];
rz(0.46485999) q[3];
sx q[3];
rz(-1.2010937) q[3];
sx q[3];
rz(0.71729898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.40806017) q[0];
sx q[0];
rz(-1.6049186) q[0];
sx q[0];
rz(-1.4419755) q[0];
rz(1.6170665) q[1];
sx q[1];
rz(-2.1433612) q[1];
sx q[1];
rz(-2.8467766) q[1];
rz(2.4649302) q[2];
sx q[2];
rz(-2.7445715) q[2];
sx q[2];
rz(-2.4286191) q[2];
rz(-0.22188998) q[3];
sx q[3];
rz(-1.1380914) q[3];
sx q[3];
rz(-1.9590845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
