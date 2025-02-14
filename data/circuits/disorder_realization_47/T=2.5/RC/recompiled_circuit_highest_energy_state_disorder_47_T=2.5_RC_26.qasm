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
rz(0.53606755) q[0];
sx q[0];
rz(-1.9789088) q[0];
sx q[0];
rz(2.7106078) q[0];
rz(0.35297901) q[1];
sx q[1];
rz(-1.9184435) q[1];
sx q[1];
rz(-2.7117742) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34944158) q[0];
sx q[0];
rz(-2.6420253) q[0];
sx q[0];
rz(-0.18886213) q[0];
rz(-pi) q[1];
rz(-1.6055029) q[2];
sx q[2];
rz(-2.7089951) q[2];
sx q[2];
rz(-0.92080599) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65548518) q[1];
sx q[1];
rz(-1.361651) q[1];
sx q[1];
rz(-0.36693348) q[1];
x q[2];
rz(1.9517035) q[3];
sx q[3];
rz(-1.105068) q[3];
sx q[3];
rz(-1.7349752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.063804403) q[2];
sx q[2];
rz(-0.52738515) q[2];
sx q[2];
rz(-1.0112313) q[2];
rz(-2.1003335) q[3];
sx q[3];
rz(-2.6285089) q[3];
sx q[3];
rz(-2.7019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1083199) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(-2.394115) q[0];
rz(2.8476818) q[1];
sx q[1];
rz(-2.0103318) q[1];
sx q[1];
rz(-0.25516137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0145484) q[0];
sx q[0];
rz(-0.60227312) q[0];
sx q[0];
rz(1.8667579) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8361144) q[2];
sx q[2];
rz(-1.7490088) q[2];
sx q[2];
rz(-2.5199948) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.98202151) q[1];
sx q[1];
rz(-2.3121093) q[1];
sx q[1];
rz(0.94653603) q[1];
x q[2];
rz(-0.18682602) q[3];
sx q[3];
rz(-2.1869724) q[3];
sx q[3];
rz(-2.0957859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1110288) q[2];
sx q[2];
rz(-0.054628987) q[2];
sx q[2];
rz(0.89152208) q[2];
rz(0.24957481) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(0.3961302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7340649) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(0.044128142) q[0];
rz(1.2491501) q[1];
sx q[1];
rz(-0.42611486) q[1];
sx q[1];
rz(0.485802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0562389) q[0];
sx q[0];
rz(-1.590343) q[0];
sx q[0];
rz(-2.9162558) q[0];
rz(2.2052231) q[2];
sx q[2];
rz(-2.9108725) q[2];
sx q[2];
rz(1.3619193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0016077) q[1];
sx q[1];
rz(-2.1473653) q[1];
sx q[1];
rz(-0.92553161) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1385192) q[3];
sx q[3];
rz(-1.2069993) q[3];
sx q[3];
rz(1.0897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8140063) q[2];
sx q[2];
rz(-1.0884716) q[2];
sx q[2];
rz(2.4468454) q[2];
rz(0.57573777) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(1.4076788) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95076743) q[0];
sx q[0];
rz(-0.29470834) q[0];
sx q[0];
rz(2.8488391) q[0];
rz(-0.5689019) q[1];
sx q[1];
rz(-2.3141373) q[1];
sx q[1];
rz(2.2946766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0305188) q[0];
sx q[0];
rz(-1.9785709) q[0];
sx q[0];
rz(0.084534377) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2402873) q[2];
sx q[2];
rz(-1.5141103) q[2];
sx q[2];
rz(-1.1336807) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49063045) q[1];
sx q[1];
rz(-1.3329026) q[1];
sx q[1];
rz(1.3722653) q[1];
rz(-1.3770695) q[3];
sx q[3];
rz(-1.5045369) q[3];
sx q[3];
rz(-0.80519262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7920821) q[2];
sx q[2];
rz(-1.258472) q[2];
sx q[2];
rz(-2.9328031) q[2];
rz(0.71792349) q[3];
sx q[3];
rz(-0.28306285) q[3];
sx q[3];
rz(0.94055241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0402886) q[0];
sx q[0];
rz(-0.52241075) q[0];
sx q[0];
rz(2.8029602) q[0];
rz(-0.70308095) q[1];
sx q[1];
rz(-2.4026726) q[1];
sx q[1];
rz(2.3032761) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0378483) q[0];
sx q[0];
rz(-2.6347343) q[0];
sx q[0];
rz(0.54308191) q[0];
rz(2.9432326) q[2];
sx q[2];
rz(-0.83399978) q[2];
sx q[2];
rz(-1.552358) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.014480249) q[1];
sx q[1];
rz(-2.4408484) q[1];
sx q[1];
rz(-1.6225918) q[1];
rz(-pi) q[2];
rz(1.1456212) q[3];
sx q[3];
rz(-2.1347649) q[3];
sx q[3];
rz(2.8759046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68359739) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(1.0028769) q[2];
rz(0.14497997) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(3.0785479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1231287) q[0];
sx q[0];
rz(-2.5557684) q[0];
sx q[0];
rz(2.1783094) q[0];
rz(-0.18114289) q[1];
sx q[1];
rz(-2.2434442) q[1];
sx q[1];
rz(-1.2394261) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9267777) q[0];
sx q[0];
rz(-2.177085) q[0];
sx q[0];
rz(-2.0303594) q[0];
rz(-pi) q[1];
rz(-2.6341637) q[2];
sx q[2];
rz(-0.52971887) q[2];
sx q[2];
rz(1.5629753) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9037348) q[1];
sx q[1];
rz(-1.9457726) q[1];
sx q[1];
rz(1.1695678) q[1];
rz(-1.419098) q[3];
sx q[3];
rz(-1.6983508) q[3];
sx q[3];
rz(-2.2383245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88508254) q[2];
sx q[2];
rz(-0.91826597) q[2];
sx q[2];
rz(-0.93290848) q[2];
rz(-1.4126623) q[3];
sx q[3];
rz(-1.4224854) q[3];
sx q[3];
rz(0.2093813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949718) q[0];
sx q[0];
rz(-0.3681204) q[0];
sx q[0];
rz(2.3387261) q[0];
rz(-2.39957) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(0.41025695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091358034) q[0];
sx q[0];
rz(-2.0977533) q[0];
sx q[0];
rz(1.9684674) q[0];
rz(3.018961) q[2];
sx q[2];
rz(-2.3735614) q[2];
sx q[2];
rz(2.1036069) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.81933016) q[1];
sx q[1];
rz(-0.17824379) q[1];
sx q[1];
rz(-2.001754) q[1];
x q[2];
rz(-2.0461844) q[3];
sx q[3];
rz(-0.96964004) q[3];
sx q[3];
rz(-0.22829311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2731169) q[2];
sx q[2];
rz(-1.5831524) q[2];
sx q[2];
rz(0.22129076) q[2];
rz(0.41915974) q[3];
sx q[3];
rz(-0.48210382) q[3];
sx q[3];
rz(-0.47590772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23685037) q[0];
sx q[0];
rz(-2.1433647) q[0];
sx q[0];
rz(1.6118443) q[0];
rz(-1.3329685) q[1];
sx q[1];
rz(-2.1905441) q[1];
sx q[1];
rz(1.453368) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4299592) q[0];
sx q[0];
rz(-0.89444133) q[0];
sx q[0];
rz(2.3948689) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4508993) q[2];
sx q[2];
rz(-1.3812307) q[2];
sx q[2];
rz(-1.3055064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9295974) q[1];
sx q[1];
rz(-2.4182165) q[1];
sx q[1];
rz(0.22774793) q[1];
rz(-1.0681719) q[3];
sx q[3];
rz(-0.31538559) q[3];
sx q[3];
rz(-2.634545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2754485) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(-0.010738372) q[2];
rz(2.8209316) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(0.58967725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43309942) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(-2.6668715) q[0];
rz(2.8482598) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(-1.6620103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728481) q[0];
sx q[0];
rz(-1.6075875) q[0];
sx q[0];
rz(-1.5307013) q[0];
rz(-pi) q[1];
rz(1.8117254) q[2];
sx q[2];
rz(-2.2325667) q[2];
sx q[2];
rz(-0.84217254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2728647) q[1];
sx q[1];
rz(-2.0365148) q[1];
sx q[1];
rz(1.3122504) q[1];
rz(2.1988499) q[3];
sx q[3];
rz(-0.30170479) q[3];
sx q[3];
rz(-2.7127271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.097229615) q[2];
sx q[2];
rz(-1.8624511) q[2];
sx q[2];
rz(-2.7105159) q[2];
rz(-2.9512682) q[3];
sx q[3];
rz(-0.35076916) q[3];
sx q[3];
rz(-2.2569807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1650319) q[0];
sx q[0];
rz(-0.29104069) q[0];
sx q[0];
rz(0.2188368) q[0];
rz(-0.95895514) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(1.3921907) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9746124) q[0];
sx q[0];
rz(-1.7638399) q[0];
sx q[0];
rz(1.6231322) q[0];
x q[1];
rz(2.8270589) q[2];
sx q[2];
rz(-2.7022903) q[2];
sx q[2];
rz(2.2404935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4569623) q[1];
sx q[1];
rz(-2.4330288) q[1];
sx q[1];
rz(2.674391) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95646071) q[3];
sx q[3];
rz(-1.9298565) q[3];
sx q[3];
rz(0.39310716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86768156) q[2];
sx q[2];
rz(-2.3751004) q[2];
sx q[2];
rz(-1.8217746) q[2];
rz(-0.53226081) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(-1.5490612) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5821447) q[0];
sx q[0];
rz(-1.5912709) q[0];
sx q[0];
rz(-2.7406319) q[0];
rz(0.22790146) q[1];
sx q[1];
rz(-2.0961998) q[1];
sx q[1];
rz(-1.6963522) q[1];
rz(-2.7264222) q[2];
sx q[2];
rz(-1.3317458) q[2];
sx q[2];
rz(1.7785704) q[2];
rz(1.8151954) q[3];
sx q[3];
rz(-2.5483589) q[3];
sx q[3];
rz(0.15729558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
