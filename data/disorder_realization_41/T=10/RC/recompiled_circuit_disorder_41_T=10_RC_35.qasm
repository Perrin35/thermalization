OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(-2.5897265) q[0];
sx q[0];
rz(-0.021835672) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(-0.2149166) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37786814) q[0];
sx q[0];
rz(-1.8452497) q[0];
sx q[0];
rz(2.8732357) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8432437) q[2];
sx q[2];
rz(-2.6539408) q[2];
sx q[2];
rz(-1.2717441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.54675198) q[1];
sx q[1];
rz(-2.0205803) q[1];
sx q[1];
rz(2.2307598) q[1];
rz(-pi) q[2];
rz(-2.4258852) q[3];
sx q[3];
rz(-1.189609) q[3];
sx q[3];
rz(0.79431278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(1.365186) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(2.2136097) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(2.2448418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0052764) q[0];
sx q[0];
rz(-1.3747842) q[0];
sx q[0];
rz(-0.85640237) q[0];
rz(2.9727544) q[2];
sx q[2];
rz(-1.4988006) q[2];
sx q[2];
rz(2.8746586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5337199) q[1];
sx q[1];
rz(-2.3405511) q[1];
sx q[1];
rz(-0.17287066) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3223022) q[3];
sx q[3];
rz(-1.2580401) q[3];
sx q[3];
rz(-2.6889338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74137694) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(2.1133912) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-0.43513402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9520156) q[0];
sx q[0];
rz(-0.36882419) q[0];
sx q[0];
rz(0.30216218) q[0];
rz(-0.52559678) q[2];
sx q[2];
rz(-0.28738775) q[2];
sx q[2];
rz(1.4488066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55089009) q[1];
sx q[1];
rz(-1.0895551) q[1];
sx q[1];
rz(2.9571556) q[1];
rz(2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(0.30291525) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(2.2241425) q[0];
rz(-0.67287412) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(0.26487574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0565856) q[0];
sx q[0];
rz(-1.6008953) q[0];
sx q[0];
rz(-2.4823275) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8633966) q[2];
sx q[2];
rz(-1.2831266) q[2];
sx q[2];
rz(-2.1759335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5323822) q[1];
sx q[1];
rz(-0.86984837) q[1];
sx q[1];
rz(2.1125395) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75986741) q[3];
sx q[3];
rz(-2.2025975) q[3];
sx q[3];
rz(2.2894273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(-1.5650361) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(-1.0331253) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(2.1889401) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8688696) q[0];
sx q[0];
rz(-0.55991828) q[0];
sx q[0];
rz(2.9193004) q[0];
rz(2.4779768) q[2];
sx q[2];
rz(-0.51082078) q[2];
sx q[2];
rz(2.7727327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2561803) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(1.0186362) q[1];
rz(-pi) q[2];
x q[2];
rz(0.080658241) q[3];
sx q[3];
rz(-0.47938743) q[3];
sx q[3];
rz(2.5731034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(0.53058132) q[2];
rz(-1.4060219) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(2.3099242) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56753165) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(-1.5166327) q[0];
rz(-1.8364871) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(2.9690202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96777746) q[0];
sx q[0];
rz(-0.41668188) q[0];
sx q[0];
rz(-1.6045051) q[0];
rz(-2.8820011) q[2];
sx q[2];
rz(-1.1278369) q[2];
sx q[2];
rz(1.358658) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.012430819) q[1];
sx q[1];
rz(-0.42814246) q[1];
sx q[1];
rz(-1.3151602) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8302912) q[3];
sx q[3];
rz(-0.65400306) q[3];
sx q[3];
rz(1.9037387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(0.78222328) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(1.8036028) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(-0.75659928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21515439) q[0];
sx q[0];
rz(-1.4190136) q[0];
sx q[0];
rz(0.093869165) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7301894) q[2];
sx q[2];
rz(-1.7559933) q[2];
sx q[2];
rz(1.411737) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8826897) q[1];
sx q[1];
rz(-1.8001302) q[1];
sx q[1];
rz(0.16112666) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6301304) q[3];
sx q[3];
rz(-1.3795128) q[3];
sx q[3];
rz(2.4997366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0044272) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(-0.19443092) q[2];
rz(-2.2284609) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(-2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(-2.1527122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522404) q[0];
sx q[0];
rz(-1.9703431) q[0];
sx q[0];
rz(-2.7826392) q[0];
rz(-3.0481911) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(1.8708558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40438548) q[1];
sx q[1];
rz(-2.7215241) q[1];
sx q[1];
rz(1.406548) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0841247) q[3];
sx q[3];
rz(-2.0351962) q[3];
sx q[3];
rz(1.7341136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6797592) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(2.3729825) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.75893629) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(-0.60920238) q[0];
rz(-3.0464879) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(2.2682155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83335919) q[0];
sx q[0];
rz(-1.5033659) q[0];
sx q[0];
rz(-1.7504577) q[0];
x q[1];
rz(1.7071502) q[2];
sx q[2];
rz(-1.7169723) q[2];
sx q[2];
rz(-1.855195) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59233353) q[1];
sx q[1];
rz(-2.1570286) q[1];
sx q[1];
rz(1.9418282) q[1];
rz(-2.14823) q[3];
sx q[3];
rz(-1.9462898) q[3];
sx q[3];
rz(1.4162228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-0.68230391) q[2];
rz(-2.7673289) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(-2.9746829) q[3];
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
rz(0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-2.1886254) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.7451161) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5293225) q[0];
sx q[0];
rz(-1.4850052) q[0];
sx q[0];
rz(-2.8180608) q[0];
rz(-pi) q[1];
rz(0.37784414) q[2];
sx q[2];
rz(-0.69381881) q[2];
sx q[2];
rz(-0.56570429) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64618387) q[1];
sx q[1];
rz(-0.42802654) q[1];
sx q[1];
rz(0.82823786) q[1];
rz(-0.20874899) q[3];
sx q[3];
rz(-2.8737846) q[3];
sx q[3];
rz(2.1684614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(2.9818025) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0970584) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(-0.13327577) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-0.940154) q[2];
sx q[2];
rz(-1.7241782) q[2];
sx q[2];
rz(0.45964514) q[2];
rz(-2.6079569) q[3];
sx q[3];
rz(-1.613637) q[3];
sx q[3];
rz(-2.2231495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
