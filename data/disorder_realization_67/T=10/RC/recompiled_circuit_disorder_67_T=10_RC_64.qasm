OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(-1.0597205) q[0];
sx q[0];
rz(0.73097316) q[0];
rz(-1.5001186) q[1];
sx q[1];
rz(-2.1067696) q[1];
sx q[1];
rz(-2.1980481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9991) q[0];
sx q[0];
rz(-1.5038135) q[0];
sx q[0];
rz(-1.5357369) q[0];
rz(0.65096345) q[2];
sx q[2];
rz(-1.5032094) q[2];
sx q[2];
rz(-2.0219321) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.01602068) q[1];
sx q[1];
rz(-1.9106094) q[1];
sx q[1];
rz(-1.1556975) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30480095) q[3];
sx q[3];
rz(-1.7497352) q[3];
sx q[3];
rz(-1.6182871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(1.250766) q[2];
rz(-1.4261774) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(-0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.132906) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(0.59894484) q[0];
rz(1.3409746) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(-0.96639955) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67993977) q[0];
sx q[0];
rz(-2.0113809) q[0];
sx q[0];
rz(-2.5011714) q[0];
rz(-pi) q[1];
rz(-1.6205377) q[2];
sx q[2];
rz(-0.94919862) q[2];
sx q[2];
rz(0.13812401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2981373) q[1];
sx q[1];
rz(-1.4569439) q[1];
sx q[1];
rz(1.650412) q[1];
rz(-0.74921272) q[3];
sx q[3];
rz(-1.2470761) q[3];
sx q[3];
rz(0.92326984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0559343) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-2.8404964) q[2];
rz(1.9484693) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10467228) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(-1.954129) q[0];
rz(-1.2359515) q[1];
sx q[1];
rz(-1.0373479) q[1];
sx q[1];
rz(1.8240066) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353969) q[0];
sx q[0];
rz(-2.4164696) q[0];
sx q[0];
rz(-2.8185185) q[0];
rz(-pi) q[1];
rz(-1.7682398) q[2];
sx q[2];
rz(-1.7806446) q[2];
sx q[2];
rz(-2.2270577) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4730075) q[1];
sx q[1];
rz(-2.7088532) q[1];
sx q[1];
rz(0.24344484) q[1];
rz(3.0242689) q[3];
sx q[3];
rz(-1.0563207) q[3];
sx q[3];
rz(2.8770212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(0.2066361) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77465039) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(2.2553717) q[0];
rz(1.0097424) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(-1.2264235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7851631) q[0];
sx q[0];
rz(-0.7298846) q[0];
sx q[0];
rz(-1.4714144) q[0];
rz(-2.1143772) q[2];
sx q[2];
rz(-2.377844) q[2];
sx q[2];
rz(2.5716512) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46982161) q[1];
sx q[1];
rz(-1.3739532) q[1];
sx q[1];
rz(-0.21398869) q[1];
rz(-pi) q[2];
rz(0.48638101) q[3];
sx q[3];
rz(-1.2380621) q[3];
sx q[3];
rz(-2.2330724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4975138) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(-2.148596) q[2];
rz(-1.3126866) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(-1.1859878) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.490373) q[1];
sx q[1];
rz(2.5591154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12594189) q[0];
sx q[0];
rz(-1.6424718) q[0];
sx q[0];
rz(-0.96836758) q[0];
rz(-pi) q[1];
x q[1];
rz(2.992606) q[2];
sx q[2];
rz(-0.89846957) q[2];
sx q[2];
rz(-0.42299262) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40618784) q[1];
sx q[1];
rz(-1.5516073) q[1];
sx q[1];
rz(2.4270942) q[1];
x q[2];
rz(1.6378239) q[3];
sx q[3];
rz(-1.1653295) q[3];
sx q[3];
rz(-1.0073347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(0.081710903) q[2];
rz(0.47406667) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(-1.4985532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24494568) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(3.1337877) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(2.0369464) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35690755) q[0];
sx q[0];
rz(-2.4570358) q[0];
sx q[0];
rz(0.25951578) q[0];
rz(-2.5014624) q[2];
sx q[2];
rz(-1.374561) q[2];
sx q[2];
rz(0.093984691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.066597477) q[1];
sx q[1];
rz(-1.0269594) q[1];
sx q[1];
rz(-0.72504136) q[1];
x q[2];
rz(-1.1052746) q[3];
sx q[3];
rz(-1.7200617) q[3];
sx q[3];
rz(2.535459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8905939) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(2.5863623) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5531439) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(0.24208367) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(-2.0297091) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2568946) q[0];
sx q[0];
rz(-1.1567133) q[0];
sx q[0];
rz(0.82722442) q[0];
rz(1.7467473) q[2];
sx q[2];
rz(-3.0627652) q[2];
sx q[2];
rz(-1.4836756) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5824077) q[1];
sx q[1];
rz(-1.8473986) q[1];
sx q[1];
rz(0.72699593) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8133409) q[3];
sx q[3];
rz(-1.8040856) q[3];
sx q[3];
rz(2.1094028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3322488) q[2];
sx q[2];
rz(-0.76433864) q[2];
sx q[2];
rz(2.3279482) q[2];
rz(1.7371197) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(0.61541921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5161045) q[0];
sx q[0];
rz(-1.3502716) q[0];
sx q[0];
rz(1.425449) q[0];
rz(-1.6199934) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-0.63751784) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.92256) q[0];
sx q[0];
rz(-2.1171283) q[0];
sx q[0];
rz(0.68740293) q[0];
rz(-pi) q[1];
rz(-2.7426863) q[2];
sx q[2];
rz(-2.3762694) q[2];
sx q[2];
rz(-1.6825636) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.70367614) q[1];
sx q[1];
rz(-2.1494467) q[1];
sx q[1];
rz(1.089536) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57070891) q[3];
sx q[3];
rz(-1.686704) q[3];
sx q[3];
rz(-1.4012208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1311538) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(-2.0054224) q[2];
rz(1.6561967) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(-1.2095399) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5159601) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(-0.67614722) q[0];
rz(2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(1.0151781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7069106) q[0];
sx q[0];
rz(-2.7043531) q[0];
sx q[0];
rz(0.31243639) q[0];
rz(1.4170309) q[2];
sx q[2];
rz(-2.2178938) q[2];
sx q[2];
rz(0.53714067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2296914) q[1];
sx q[1];
rz(-1.9209314) q[1];
sx q[1];
rz(1.196911) q[1];
rz(-pi) q[2];
rz(-2.3750651) q[3];
sx q[3];
rz(-1.3228647) q[3];
sx q[3];
rz(2.8299696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(2.1742415) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5230781) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.4046232) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.343095) q[0];
sx q[0];
rz(-1.6999177) q[0];
sx q[0];
rz(-1.9341747) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8374412) q[2];
sx q[2];
rz(-1.6870105) q[2];
sx q[2];
rz(-0.2955557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9060668) q[1];
sx q[1];
rz(-2.0350254) q[1];
sx q[1];
rz(2.3266351) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9910562) q[3];
sx q[3];
rz(-1.5621462) q[3];
sx q[3];
rz(2.0899525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.848032) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(-2.8354697) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(-1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.33481471) q[0];
sx q[0];
rz(-2.0756742) q[0];
sx q[0];
rz(-2.6609127) q[0];
rz(-1.9819992) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(-1.8691312) q[2];
sx q[2];
rz(-1.3186426) q[2];
sx q[2];
rz(-1.1124055) q[2];
rz(2.0648099) q[3];
sx q[3];
rz(-1.2772588) q[3];
sx q[3];
rz(0.96989934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
