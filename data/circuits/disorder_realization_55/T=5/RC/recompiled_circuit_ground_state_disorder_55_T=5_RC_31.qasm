OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2321229) q[0];
sx q[0];
rz(2.1535518) q[0];
sx q[0];
rz(9.0721985) q[0];
rz(2.4352788) q[1];
sx q[1];
rz(-2.3568454) q[1];
sx q[1];
rz(3.1142601) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407878) q[0];
sx q[0];
rz(-0.39310021) q[0];
sx q[0];
rz(2.604305) q[0];
x q[1];
rz(0.17345239) q[2];
sx q[2];
rz(-1.7043742) q[2];
sx q[2];
rz(1.1277778) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24522745) q[1];
sx q[1];
rz(-1.7385529) q[1];
sx q[1];
rz(-1.6025402) q[1];
rz(-0.54630791) q[3];
sx q[3];
rz(-0.27537333) q[3];
sx q[3];
rz(0.34378036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0722384) q[2];
sx q[2];
rz(-1.8202929) q[2];
sx q[2];
rz(-1.5864774) q[2];
rz(0.72748264) q[3];
sx q[3];
rz(-0.95557094) q[3];
sx q[3];
rz(2.6684842) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.469406) q[0];
sx q[0];
rz(-1.7962026) q[0];
sx q[0];
rz(0.94959062) q[0];
rz(0.32497111) q[1];
sx q[1];
rz(-1.3409216) q[1];
sx q[1];
rz(2.6281338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9531813) q[0];
sx q[0];
rz(-1.1674178) q[0];
sx q[0];
rz(-0.90103407) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32294257) q[2];
sx q[2];
rz(-2.2720798) q[2];
sx q[2];
rz(2.8632134) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1850441) q[1];
sx q[1];
rz(-2.0868851) q[1];
sx q[1];
rz(0.98770492) q[1];
x q[2];
rz(-0.49118941) q[3];
sx q[3];
rz(-1.0460539) q[3];
sx q[3];
rz(-2.3150746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.010633858) q[2];
sx q[2];
rz(-1.8124688) q[2];
sx q[2];
rz(-2.0531674) q[2];
rz(-2.4736577) q[3];
sx q[3];
rz(-0.1440983) q[3];
sx q[3];
rz(-1.4151423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.3208991) q[0];
sx q[0];
rz(-2.2405393) q[0];
sx q[0];
rz(-1.2574842) q[0];
rz(-3.056774) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(0.65748293) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5710058) q[0];
sx q[0];
rz(-2.0181542) q[0];
sx q[0];
rz(1.4259095) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74983238) q[2];
sx q[2];
rz(-0.74374357) q[2];
sx q[2];
rz(2.2735571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8264933) q[1];
sx q[1];
rz(-1.4001453) q[1];
sx q[1];
rz(2.6572589) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30835521) q[3];
sx q[3];
rz(-0.93290795) q[3];
sx q[3];
rz(-2.9264192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.377044) q[2];
sx q[2];
rz(-2.720764) q[2];
sx q[2];
rz(1.8910889) q[2];
rz(1.6999543) q[3];
sx q[3];
rz(-1.7504033) q[3];
sx q[3];
rz(0.82956782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4568951) q[0];
sx q[0];
rz(-0.95624113) q[0];
sx q[0];
rz(2.5392505) q[0];
rz(0.96386987) q[1];
sx q[1];
rz(-2.3606221) q[1];
sx q[1];
rz(0.22183713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9070061) q[0];
sx q[0];
rz(-0.013628634) q[0];
sx q[0];
rz(0.6032633) q[0];
x q[1];
rz(-2.8922563) q[2];
sx q[2];
rz(-2.0866924) q[2];
sx q[2];
rz(-0.14125401) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.79961046) q[1];
sx q[1];
rz(-2.1253808) q[1];
sx q[1];
rz(1.0971402) q[1];
rz(-pi) q[2];
rz(-2.6915914) q[3];
sx q[3];
rz(-1.5061989) q[3];
sx q[3];
rz(0.76417506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43231371) q[2];
sx q[2];
rz(-1.9807434) q[2];
sx q[2];
rz(-0.088689001) q[2];
rz(-0.49897075) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(-3.0626845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0093507) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(0.90231878) q[0];
rz(0.29014507) q[1];
sx q[1];
rz(-2.2369308) q[1];
sx q[1];
rz(-1.1660928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48448823) q[0];
sx q[0];
rz(-2.1323626) q[0];
sx q[0];
rz(-1.7620371) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71641123) q[2];
sx q[2];
rz(-0.36135095) q[2];
sx q[2];
rz(0.30308576) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.13671) q[1];
sx q[1];
rz(-0.45082475) q[1];
sx q[1];
rz(0.87597998) q[1];
x q[2];
rz(-2.8899447) q[3];
sx q[3];
rz(-1.5960428) q[3];
sx q[3];
rz(1.3803052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0151998) q[2];
sx q[2];
rz(-2.3455399) q[2];
sx q[2];
rz(-1.2882721) q[2];
rz(-2.1081693) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(-1.4748632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7159336) q[0];
sx q[0];
rz(-0.10361828) q[0];
sx q[0];
rz(2.9820251) q[0];
rz(2.5345934) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(1.0608231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95283857) q[0];
sx q[0];
rz(-2.5079281) q[0];
sx q[0];
rz(1.2418644) q[0];
rz(-2.873704) q[2];
sx q[2];
rz(-1.3527591) q[2];
sx q[2];
rz(-1.5341369) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2810106) q[1];
sx q[1];
rz(-1.7063024) q[1];
sx q[1];
rz(2.1675088) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6688329) q[3];
sx q[3];
rz(-0.99409249) q[3];
sx q[3];
rz(-2.2876379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0173505) q[2];
sx q[2];
rz(-2.3257747) q[2];
sx q[2];
rz(-0.38786495) q[2];
rz(-1.3346437) q[3];
sx q[3];
rz(-0.67238656) q[3];
sx q[3];
rz(2.1147125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.89011985) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(0.89333308) q[0];
rz(2.6003301) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(-2.0792686) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5389299) q[0];
sx q[0];
rz(-1.3100123) q[0];
sx q[0];
rz(-2.3988924) q[0];
x q[1];
rz(2.3184653) q[2];
sx q[2];
rz(-1.6170653) q[2];
sx q[2];
rz(-2.0101765) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3596284) q[1];
sx q[1];
rz(-2.5853695) q[1];
sx q[1];
rz(3.0969949) q[1];
rz(-pi) q[2];
rz(-1.7797864) q[3];
sx q[3];
rz(-2.3694443) q[3];
sx q[3];
rz(1.0812372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58747753) q[2];
sx q[2];
rz(-0.97439659) q[2];
sx q[2];
rz(1.1678196) q[2];
rz(-2.3217412) q[3];
sx q[3];
rz(-2.3746115) q[3];
sx q[3];
rz(-2.4988373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.5646566) q[0];
sx q[0];
rz(-0.37207237) q[0];
sx q[0];
rz(-1.8336953) q[0];
rz(-1.4564184) q[1];
sx q[1];
rz(-1.6273727) q[1];
sx q[1];
rz(0.13052043) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3878138) q[0];
sx q[0];
rz(-2.8243833) q[0];
sx q[0];
rz(-1.3340201) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0575664) q[2];
sx q[2];
rz(-2.3064838) q[2];
sx q[2];
rz(-1.4722919) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.710142) q[1];
sx q[1];
rz(-1.6098861) q[1];
sx q[1];
rz(-2.4907656) q[1];
x q[2];
rz(-0.91754387) q[3];
sx q[3];
rz(-0.96644478) q[3];
sx q[3];
rz(1.0928327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.1859583) q[2];
sx q[2];
rz(-2.3027577) q[2];
sx q[2];
rz(-1.4107417) q[2];
rz(-0.12061067) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(1.6284778) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2738344) q[0];
sx q[0];
rz(-2.5262316) q[0];
sx q[0];
rz(0.14061418) q[0];
rz(0.73453844) q[1];
sx q[1];
rz(-1.4334375) q[1];
sx q[1];
rz(2.3908884) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1575506) q[0];
sx q[0];
rz(-2.7617402) q[0];
sx q[0];
rz(-2.4628839) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23739135) q[2];
sx q[2];
rz(-1.2331881) q[2];
sx q[2];
rz(-2.2142976) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6748811) q[1];
sx q[1];
rz(-1.3916236) q[1];
sx q[1];
rz(3.0790943) q[1];
rz(0.073551579) q[3];
sx q[3];
rz(-2.2760609) q[3];
sx q[3];
rz(1.2770231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1371896) q[2];
sx q[2];
rz(-2.4418094) q[2];
sx q[2];
rz(2.3889551) q[2];
rz(0.33958069) q[3];
sx q[3];
rz(-1.3656253) q[3];
sx q[3];
rz(-3.1201709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4362157) q[0];
sx q[0];
rz(-0.74420539) q[0];
sx q[0];
rz(2.5035653) q[0];
rz(2.4127507) q[1];
sx q[1];
rz(-1.809027) q[1];
sx q[1];
rz(2.3400838) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1164393) q[0];
sx q[0];
rz(-2.354216) q[0];
sx q[0];
rz(1.8579432) q[0];
rz(2.3036868) q[2];
sx q[2];
rz(-1.3896754) q[2];
sx q[2];
rz(1.6194026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.90097016) q[1];
sx q[1];
rz(-0.83604807) q[1];
sx q[1];
rz(0.052608629) q[1];
rz(-pi) q[2];
rz(1.784104) q[3];
sx q[3];
rz(-0.65591988) q[3];
sx q[3];
rz(-1.8260969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8375497) q[2];
sx q[2];
rz(-1.7767228) q[2];
sx q[2];
rz(1.1207646) q[2];
rz(2.4325727) q[3];
sx q[3];
rz(-0.46375436) q[3];
sx q[3];
rz(1.2254865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7365702) q[0];
sx q[0];
rz(-0.46157349) q[0];
sx q[0];
rz(2.9622958) q[0];
rz(-0.90514056) q[1];
sx q[1];
rz(-0.97733472) q[1];
sx q[1];
rz(2.6505145) q[1];
rz(0.84409406) q[2];
sx q[2];
rz(-1.8432968) q[2];
sx q[2];
rz(-2.7805614) q[2];
rz(2.8560588) q[3];
sx q[3];
rz(-2.0750506) q[3];
sx q[3];
rz(-1.688523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
