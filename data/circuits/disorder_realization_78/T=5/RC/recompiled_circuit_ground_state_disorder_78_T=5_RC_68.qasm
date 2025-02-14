OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8371589) q[0];
sx q[0];
rz(4.9871939) q[0];
sx q[0];
rz(11.468588) q[0];
rz(2.1863565) q[1];
sx q[1];
rz(-0.74164852) q[1];
sx q[1];
rz(0.15129605) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77028217) q[0];
sx q[0];
rz(-1.4454968) q[0];
sx q[0];
rz(-0.025733982) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.120235) q[2];
sx q[2];
rz(-1.6701227) q[2];
sx q[2];
rz(-1.8252357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23257593) q[1];
sx q[1];
rz(-2.1409509) q[1];
sx q[1];
rz(-2.2297165) q[1];
x q[2];
rz(2.339974) q[3];
sx q[3];
rz(-1.3676757) q[3];
sx q[3];
rz(2.5020848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1260881) q[2];
sx q[2];
rz(-1.2903004) q[2];
sx q[2];
rz(0.31910953) q[2];
rz(2.0305521) q[3];
sx q[3];
rz(-0.55912656) q[3];
sx q[3];
rz(1.3148974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1183209) q[0];
sx q[0];
rz(-0.51902223) q[0];
sx q[0];
rz(0.58854377) q[0];
rz(0.59666807) q[1];
sx q[1];
rz(-1.3304973) q[1];
sx q[1];
rz(-0.24135022) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7241192) q[0];
sx q[0];
rz(-2.3406174) q[0];
sx q[0];
rz(-2.1885314) q[0];
rz(-pi) q[1];
rz(-1.4407519) q[2];
sx q[2];
rz(-0.50419129) q[2];
sx q[2];
rz(0.90182226) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4171364) q[1];
sx q[1];
rz(-2.0057185) q[1];
sx q[1];
rz(-3.0318854) q[1];
rz(-pi) q[2];
rz(0.94542687) q[3];
sx q[3];
rz(-1.4875879) q[3];
sx q[3];
rz(2.0997467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1179463) q[2];
sx q[2];
rz(-1.6829374) q[2];
sx q[2];
rz(-3.0493951) q[2];
rz(-2.2495031) q[3];
sx q[3];
rz(-0.8725608) q[3];
sx q[3];
rz(2.0766855) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43111619) q[0];
sx q[0];
rz(-1.5065864) q[0];
sx q[0];
rz(-3.0431252) q[0];
rz(-0.59421986) q[1];
sx q[1];
rz(-2.0829945) q[1];
sx q[1];
rz(2.1509511) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0069283) q[0];
sx q[0];
rz(-1.1968839) q[0];
sx q[0];
rz(-0.093061826) q[0];
rz(-pi) q[1];
rz(-2.1192693) q[2];
sx q[2];
rz(-2.2532941) q[2];
sx q[2];
rz(-2.2816531) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16690635) q[1];
sx q[1];
rz(-2.1975127) q[1];
sx q[1];
rz(0.30345602) q[1];
rz(-1.1029937) q[3];
sx q[3];
rz(-1.3672223) q[3];
sx q[3];
rz(0.14467231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46382612) q[2];
sx q[2];
rz(-0.78654424) q[2];
sx q[2];
rz(-0.80219913) q[2];
rz(-1.5944611) q[3];
sx q[3];
rz(-2.0764669) q[3];
sx q[3];
rz(-2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7441854) q[0];
sx q[0];
rz(-0.35905251) q[0];
sx q[0];
rz(-1.3209976) q[0];
rz(1.6162704) q[1];
sx q[1];
rz(-0.95322144) q[1];
sx q[1];
rz(0.1604518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67293834) q[0];
sx q[0];
rz(-0.78862353) q[0];
sx q[0];
rz(1.9140052) q[0];
rz(0.58893369) q[2];
sx q[2];
rz(-2.9351165) q[2];
sx q[2];
rz(1.4331417) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0640434) q[1];
sx q[1];
rz(-2.4740015) q[1];
sx q[1];
rz(1.8197219) q[1];
rz(1.1536648) q[3];
sx q[3];
rz(-2.0685745) q[3];
sx q[3];
rz(-2.1912632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3204331) q[2];
sx q[2];
rz(-1.4904138) q[2];
sx q[2];
rz(-2.5725345) q[2];
rz(-1.9197561) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(3.0088185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4959167) q[0];
sx q[0];
rz(-0.78302947) q[0];
sx q[0];
rz(0.83086479) q[0];
rz(-1.2105385) q[1];
sx q[1];
rz(-1.6485063) q[1];
sx q[1];
rz(-0.99162203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094153397) q[0];
sx q[0];
rz(-1.2625361) q[0];
sx q[0];
rz(0.83931132) q[0];
rz(-1.2328202) q[2];
sx q[2];
rz(-2.4261279) q[2];
sx q[2];
rz(-1.6140661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8960921) q[1];
sx q[1];
rz(-0.84907167) q[1];
sx q[1];
rz(0.20770276) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5414821) q[3];
sx q[3];
rz(-0.94621822) q[3];
sx q[3];
rz(1.2219714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7603989) q[2];
sx q[2];
rz(-0.92279592) q[2];
sx q[2];
rz(-2.7823616) q[2];
rz(2.4257816) q[3];
sx q[3];
rz(-2.3575213) q[3];
sx q[3];
rz(1.1904967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0773709) q[0];
sx q[0];
rz(-1.2718028) q[0];
sx q[0];
rz(2.6203058) q[0];
rz(2.298666) q[1];
sx q[1];
rz(-1.9631674) q[1];
sx q[1];
rz(3.0311323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9006289) q[0];
sx q[0];
rz(-0.74375737) q[0];
sx q[0];
rz(-2.7636011) q[0];
rz(-2.213352) q[2];
sx q[2];
rz(-1.5470328) q[2];
sx q[2];
rz(1.0500963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46347324) q[1];
sx q[1];
rz(-1.7686525) q[1];
sx q[1];
rz(1.0956647) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0024226) q[3];
sx q[3];
rz(-2.5768498) q[3];
sx q[3];
rz(3.125001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20737401) q[2];
sx q[2];
rz(-1.5418345) q[2];
sx q[2];
rz(-2.1477487) q[2];
rz(0.55073109) q[3];
sx q[3];
rz(-2.6271074) q[3];
sx q[3];
rz(-1.6186835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.9024502) q[0];
sx q[0];
rz(-0.37343326) q[0];
sx q[0];
rz(-0.19110876) q[0];
rz(2.7725819) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(-0.63327995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8610982) q[0];
sx q[0];
rz(-1.3320574) q[0];
sx q[0];
rz(-1.6709836) q[0];
rz(0.72004135) q[2];
sx q[2];
rz(-1.9332464) q[2];
sx q[2];
rz(1.8099648) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7821473) q[1];
sx q[1];
rz(-1.7759062) q[1];
sx q[1];
rz(-2.9444749) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4426887) q[3];
sx q[3];
rz(-1.4443384) q[3];
sx q[3];
rz(1.5497643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9621027) q[2];
sx q[2];
rz(-1.3733764) q[2];
sx q[2];
rz(-2.7607259) q[2];
rz(1.7535836) q[3];
sx q[3];
rz(-1.0264779) q[3];
sx q[3];
rz(-1.7109722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66985828) q[0];
sx q[0];
rz(-2.7405881) q[0];
sx q[0];
rz(-1.5420472) q[0];
rz(1.3767287) q[1];
sx q[1];
rz(-1.7574666) q[1];
sx q[1];
rz(-0.9309887) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.72351) q[0];
sx q[0];
rz(-0.93659725) q[0];
sx q[0];
rz(-1.67976) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8225708) q[2];
sx q[2];
rz(-2.012017) q[2];
sx q[2];
rz(-0.069442858) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38057391) q[1];
sx q[1];
rz(-2.3133458) q[1];
sx q[1];
rz(3.1395802) q[1];
x q[2];
rz(-2.2119207) q[3];
sx q[3];
rz(-0.76477178) q[3];
sx q[3];
rz(-1.1906033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5652183) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(-3.0726688) q[2];
rz(1.8779514) q[3];
sx q[3];
rz(-2.7694323) q[3];
sx q[3];
rz(-2.4431958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3987592) q[0];
sx q[0];
rz(-2.1837406) q[0];
sx q[0];
rz(-3.0897019) q[0];
rz(0.17008153) q[1];
sx q[1];
rz(-2.4982128) q[1];
sx q[1];
rz(-0.35194078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1387685) q[0];
sx q[0];
rz(-1.616298) q[0];
sx q[0];
rz(2.7753434) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62656709) q[2];
sx q[2];
rz(-2.0800903) q[2];
sx q[2];
rz(1.9208822) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.62427795) q[1];
sx q[1];
rz(-1.0464962) q[1];
sx q[1];
rz(0.45714) q[1];
x q[2];
rz(-1.264335) q[3];
sx q[3];
rz(-1.500794) q[3];
sx q[3];
rz(0.7037735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45741442) q[2];
sx q[2];
rz(-2.4922721) q[2];
sx q[2];
rz(-2.7049098) q[2];
rz(-0.39143482) q[3];
sx q[3];
rz(-1.6792363) q[3];
sx q[3];
rz(2.2552538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201685) q[0];
sx q[0];
rz(-2.4346209) q[0];
sx q[0];
rz(-0.11908764) q[0];
rz(1.2991615) q[1];
sx q[1];
rz(-1.7487339) q[1];
sx q[1];
rz(1.3577168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8004476) q[0];
sx q[0];
rz(-0.94609944) q[0];
sx q[0];
rz(-0.50826061) q[0];
x q[1];
rz(2.7804271) q[2];
sx q[2];
rz(-2.2884011) q[2];
sx q[2];
rz(1.4360365) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5636041) q[1];
sx q[1];
rz(-1.9050026) q[1];
sx q[1];
rz(-0.6329221) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81858738) q[3];
sx q[3];
rz(-1.162611) q[3];
sx q[3];
rz(-3.0455923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9327717) q[2];
sx q[2];
rz(-1.5522771) q[2];
sx q[2];
rz(2.5271752) q[2];
rz(1.6299853) q[3];
sx q[3];
rz(-2.4698518) q[3];
sx q[3];
rz(0.083812788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3096302) q[0];
sx q[0];
rz(-0.93292581) q[0];
sx q[0];
rz(0.25136872) q[0];
rz(-2.108719) q[1];
sx q[1];
rz(-1.1943457) q[1];
sx q[1];
rz(1.7051382) q[1];
rz(-0.9821427) q[2];
sx q[2];
rz(-0.12599421) q[2];
sx q[2];
rz(-0.85859184) q[2];
rz(2.4182416) q[3];
sx q[3];
rz(-0.94905973) q[3];
sx q[3];
rz(2.2524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
