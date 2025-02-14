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
rz(0.41843721) q[0];
sx q[0];
rz(-0.96324459) q[0];
sx q[0];
rz(0.20382398) q[0];
rz(1.0526429) q[1];
sx q[1];
rz(3.9349603) q[1];
sx q[1];
rz(10.959672) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14102916) q[0];
sx q[0];
rz(-1.3357497) q[0];
sx q[0];
rz(-1.1147333) q[0];
rz(-pi) q[1];
rz(-0.84759961) q[2];
sx q[2];
rz(-1.3061451) q[2];
sx q[2];
rz(-2.6516595) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4272449) q[1];
sx q[1];
rz(-0.74810076) q[1];
sx q[1];
rz(1.0393049) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9382557) q[3];
sx q[3];
rz(-1.4383738) q[3];
sx q[3];
rz(-1.2156957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4172998) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(2.144045) q[2];
rz(-0.7558465) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(-2.1006179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7928829) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(-1.8541699) q[0];
rz(-2.6584794) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(1.2189254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12139509) q[0];
sx q[0];
rz(-1.7292062) q[0];
sx q[0];
rz(-1.3729457) q[0];
x q[1];
rz(2.4697815) q[2];
sx q[2];
rz(-1.3992568) q[2];
sx q[2];
rz(1.1900657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9049553) q[1];
sx q[1];
rz(-1.7153011) q[1];
sx q[1];
rz(-1.852024) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7876625) q[3];
sx q[3];
rz(-0.90995379) q[3];
sx q[3];
rz(-2.5315745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7634742) q[2];
sx q[2];
rz(-3.0749622) q[2];
sx q[2];
rz(-2.5197869) q[2];
rz(0.63831896) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8353552) q[0];
sx q[0];
rz(-0.64202809) q[0];
sx q[0];
rz(-2.9252692) q[0];
rz(-0.59858876) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(-2.6187706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.110958) q[0];
sx q[0];
rz(-1.2485663) q[0];
sx q[0];
rz(2.8186174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4141358) q[2];
sx q[2];
rz(-1.1839424) q[2];
sx q[2];
rz(-1.1413107) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.712341) q[1];
sx q[1];
rz(-2.0852226) q[1];
sx q[1];
rz(1.3144668) q[1];
rz(0.30662068) q[3];
sx q[3];
rz(-1.8869747) q[3];
sx q[3];
rz(0.090360377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.181695) q[2];
sx q[2];
rz(-1.8671702) q[2];
sx q[2];
rz(1.6240906) q[2];
rz(2.6767139) q[3];
sx q[3];
rz(-2.2113694) q[3];
sx q[3];
rz(1.6200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14438039) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(1.5390747) q[0];
rz(2.1082361) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(-1.296952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3281477) q[0];
sx q[0];
rz(-1.362365) q[0];
sx q[0];
rz(2.3013902) q[0];
rz(1.4960852) q[2];
sx q[2];
rz(-2.4444408) q[2];
sx q[2];
rz(2.8139204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32758157) q[1];
sx q[1];
rz(-1.1594698) q[1];
sx q[1];
rz(-0.41366215) q[1];
rz(-pi) q[2];
rz(0.4850895) q[3];
sx q[3];
rz(-0.61044932) q[3];
sx q[3];
rz(1.7318673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79232717) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(2.0959334) q[2];
rz(0.21444923) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8056718) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(2.6222498) q[0];
rz(0.94201159) q[1];
sx q[1];
rz(-0.28128925) q[1];
sx q[1];
rz(-1.5325783) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2359152) q[0];
sx q[0];
rz(-0.92936838) q[0];
sx q[0];
rz(-1.1113313) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4283044) q[2];
sx q[2];
rz(-1.7497369) q[2];
sx q[2];
rz(-1.475309) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11621257) q[1];
sx q[1];
rz(-1.0702225) q[1];
sx q[1];
rz(-0.0059555014) q[1];
rz(-0.71096731) q[3];
sx q[3];
rz(-0.14440726) q[3];
sx q[3];
rz(-0.84916964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64451009) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(-0.28437781) q[2];
rz(-2.583368) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(-0.24793454) q[3];
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
rz(0.068950653) q[0];
sx q[0];
rz(-2.729029) q[0];
sx q[0];
rz(-0.64055842) q[0];
rz(1.5156281) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(0.21387771) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2757077) q[0];
sx q[0];
rz(-1.2684221) q[0];
sx q[0];
rz(-2.9024966) q[0];
rz(0.40753813) q[2];
sx q[2];
rz(-0.78754163) q[2];
sx q[2];
rz(2.119273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8941514) q[1];
sx q[1];
rz(-0.94235984) q[1];
sx q[1];
rz(-0.085606023) q[1];
rz(-pi) q[2];
rz(1.050726) q[3];
sx q[3];
rz(-1.5087391) q[3];
sx q[3];
rz(2.7620535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3512257) q[2];
sx q[2];
rz(-2.6595317) q[2];
sx q[2];
rz(-2.9638929) q[2];
rz(-1.3972345) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(-2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099667065) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(-1.1055111) q[0];
rz(-0.28327495) q[1];
sx q[1];
rz(-0.53422821) q[1];
sx q[1];
rz(-1.6800605) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7992226) q[0];
sx q[0];
rz(-2.2307751) q[0];
sx q[0];
rz(-2.2837385) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2143557) q[2];
sx q[2];
rz(-0.74802784) q[2];
sx q[2];
rz(3.1035556) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0431598) q[1];
sx q[1];
rz(-0.57505703) q[1];
sx q[1];
rz(-0.21456031) q[1];
rz(1.2456149) q[3];
sx q[3];
rz(-1.0187314) q[3];
sx q[3];
rz(1.0528885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0854411) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(-0.20480569) q[2];
rz(1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-0.36627305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751223) q[0];
sx q[0];
rz(-0.46638745) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(-1.9974476) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9947203) q[0];
sx q[0];
rz(-1.084274) q[0];
sx q[0];
rz(-1.542132) q[0];
x q[1];
rz(1.8587684) q[2];
sx q[2];
rz(-0.33703732) q[2];
sx q[2];
rz(2.2971643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.856196) q[1];
sx q[1];
rz(-0.90803972) q[1];
sx q[1];
rz(1.0876571) q[1];
rz(1.2253292) q[3];
sx q[3];
rz(-1.8356712) q[3];
sx q[3];
rz(0.027907413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3860151) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(-1.8899274) q[2];
rz(-1.7898611) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(-2.1613817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0018472483) q[0];
sx q[0];
rz(-0.63848764) q[0];
sx q[0];
rz(-1.2563323) q[0];
rz(0.095257692) q[1];
sx q[1];
rz(-2.2525411) q[1];
sx q[1];
rz(-2.6108066) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4149218) q[0];
sx q[0];
rz(-1.7309395) q[0];
sx q[0];
rz(1.8216351) q[0];
x q[1];
rz(0.31802788) q[2];
sx q[2];
rz(-1.6075168) q[2];
sx q[2];
rz(1.5098234) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7698313) q[1];
sx q[1];
rz(-0.27250817) q[1];
sx q[1];
rz(-2.647577) q[1];
rz(-pi) q[2];
rz(-1.68566) q[3];
sx q[3];
rz(-1.126976) q[3];
sx q[3];
rz(0.97568363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8119767) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(0.54083332) q[2];
rz(-2.3234308) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(-2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.624991) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(-0.3987819) q[0];
rz(-1.3840236) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(0.049364518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2731664) q[0];
sx q[0];
rz(-3.0169562) q[0];
sx q[0];
rz(-1.0602555) q[0];
rz(1.3081506) q[2];
sx q[2];
rz(-1.4672973) q[2];
sx q[2];
rz(-1.7067133) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2422707) q[1];
sx q[1];
rz(-0.99192109) q[1];
sx q[1];
rz(3.088356) q[1];
rz(-pi) q[2];
rz(0.27487288) q[3];
sx q[3];
rz(-2.4913553) q[3];
sx q[3];
rz(0.45675983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0177239) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(-1.5105985) q[2];
rz(2.2137568) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.8265726) q[0];
sx q[0];
rz(-0.44354225) q[0];
sx q[0];
rz(-1.0916239) q[0];
rz(1.0646959) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(2.3774316) q[2];
sx q[2];
rz(-2.4780826) q[2];
sx q[2];
rz(1.6586951) q[2];
rz(1.0743027) q[3];
sx q[3];
rz(-0.60582325) q[3];
sx q[3];
rz(-0.7614991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
