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
rz(3.1177899) q[0];
sx q[0];
rz(-1.1060214) q[0];
sx q[0];
rz(-0.7769146) q[0];
rz(-1.7493526) q[1];
sx q[1];
rz(-1.8267781) q[1];
sx q[1];
rz(-2.1652752) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38589729) q[0];
sx q[0];
rz(-1.160826) q[0];
sx q[0];
rz(0.14571054) q[0];
rz(-pi) q[1];
rz(0.56371477) q[2];
sx q[2];
rz(-1.6275121) q[2];
sx q[2];
rz(-2.4617755) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7853422) q[1];
sx q[1];
rz(-1.1467198) q[1];
sx q[1];
rz(-0.29845684) q[1];
rz(-pi) q[2];
x q[2];
rz(0.032194897) q[3];
sx q[3];
rz(-2.1436084) q[3];
sx q[3];
rz(0.16530748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6355847) q[2];
sx q[2];
rz(-0.053746544) q[2];
sx q[2];
rz(-0.40679833) q[2];
rz(-0.16945101) q[3];
sx q[3];
rz(-2.6120766) q[3];
sx q[3];
rz(1.0725526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605543) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(-0.0037923092) q[0];
rz(0.077839851) q[1];
sx q[1];
rz(-0.65996116) q[1];
sx q[1];
rz(-0.30581623) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77340375) q[0];
sx q[0];
rz(-2.0819252) q[0];
sx q[0];
rz(-0.21171932) q[0];
rz(-0.75188101) q[2];
sx q[2];
rz(-1.6691748) q[2];
sx q[2];
rz(1.6461314) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9287323) q[1];
sx q[1];
rz(-1.1747735) q[1];
sx q[1];
rz(-0.81800445) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3313794) q[3];
sx q[3];
rz(-3.0533724) q[3];
sx q[3];
rz(2.719413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1514312) q[2];
sx q[2];
rz(-1.7502681) q[2];
sx q[2];
rz(0.64275098) q[2];
rz(3.0691872) q[3];
sx q[3];
rz(-2.0824771) q[3];
sx q[3];
rz(-1.7806627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533503) q[0];
sx q[0];
rz(-2.998816) q[0];
sx q[0];
rz(0.15637583) q[0];
rz(0.0414255) q[1];
sx q[1];
rz(-0.62774575) q[1];
sx q[1];
rz(-1.590439) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.56083) q[0];
sx q[0];
rz(-1.3645882) q[0];
sx q[0];
rz(-2.5872562) q[0];
rz(1.6386436) q[2];
sx q[2];
rz(-0.91582752) q[2];
sx q[2];
rz(-2.613435) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.157718) q[1];
sx q[1];
rz(-1.0441171) q[1];
sx q[1];
rz(1.7009363) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66371347) q[3];
sx q[3];
rz(-0.98582375) q[3];
sx q[3];
rz(3.0970517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.40858832) q[2];
sx q[2];
rz(-0.53885794) q[2];
sx q[2];
rz(3.1206701) q[2];
rz(-2.9544592) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(-0.11370295) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1784096) q[0];
sx q[0];
rz(-2.8091176) q[0];
sx q[0];
rz(-2.7030429) q[0];
rz(-1.5248388) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(-2.8964892) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.995718) q[0];
sx q[0];
rz(-2.2882713) q[0];
sx q[0];
rz(0.11086734) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4133873) q[2];
sx q[2];
rz(-1.206996) q[2];
sx q[2];
rz(-0.25870332) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7286543) q[1];
sx q[1];
rz(-2.1082343) q[1];
sx q[1];
rz(-2.9793903) q[1];
rz(2.80527) q[3];
sx q[3];
rz(-2.3592279) q[3];
sx q[3];
rz(2.4049644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(2.8098246) q[2];
rz(-0.48745421) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(0.99307466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29193923) q[0];
sx q[0];
rz(-1.6878457) q[0];
sx q[0];
rz(-0.77350235) q[0];
rz(-1.1812814) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(-1.7519417) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6365125) q[0];
sx q[0];
rz(-2.5861222) q[0];
sx q[0];
rz(2.4019) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2393786) q[2];
sx q[2];
rz(-2.3093975) q[2];
sx q[2];
rz(0.26022831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2323245) q[1];
sx q[1];
rz(-1.0977912) q[1];
sx q[1];
rz(1.478341) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23791194) q[3];
sx q[3];
rz(-1.7470659) q[3];
sx q[3];
rz(0.59559866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2191849) q[2];
sx q[2];
rz(-1.2097404) q[2];
sx q[2];
rz(2.6694471) q[2];
rz(-1.8426497) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(-0.81056547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.2209114) q[0];
sx q[0];
rz(-0.26316106) q[0];
sx q[0];
rz(0.26350185) q[0];
rz(-1.1031411) q[1];
sx q[1];
rz(-1.3145072) q[1];
sx q[1];
rz(-0.37364328) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140681) q[0];
sx q[0];
rz(-1.4980982) q[0];
sx q[0];
rz(-0.060427314) q[0];
rz(-1.4336072) q[2];
sx q[2];
rz(-2.409916) q[2];
sx q[2];
rz(-0.2889932) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3195575) q[1];
sx q[1];
rz(-0.6614092) q[1];
sx q[1];
rz(1.0378077) q[1];
x q[2];
rz(-0.75833894) q[3];
sx q[3];
rz(-0.74303526) q[3];
sx q[3];
rz(2.5287927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11334795) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(-2.587758) q[2];
rz(1.7438186) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(2.8288614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5572307) q[0];
sx q[0];
rz(-2.1350242) q[0];
sx q[0];
rz(-1.2297909) q[0];
rz(-0.23904414) q[1];
sx q[1];
rz(-1.5115279) q[1];
sx q[1];
rz(2.8412433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8431138) q[0];
sx q[0];
rz(-1.4550147) q[0];
sx q[0];
rz(-2.9779469) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5671674) q[2];
sx q[2];
rz(-2.3544899) q[2];
sx q[2];
rz(-1.8970823) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9508267) q[1];
sx q[1];
rz(-3.1263104) q[1];
sx q[1];
rz(-1.5453668) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85791479) q[3];
sx q[3];
rz(-1.5893717) q[3];
sx q[3];
rz(-1.9481079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.131669) q[2];
sx q[2];
rz(-1.6841623) q[2];
sx q[2];
rz(-0.24492502) q[2];
rz(2.6217672) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(2.4533217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7898665) q[0];
sx q[0];
rz(-0.34857294) q[0];
sx q[0];
rz(1.9785731) q[0];
rz(3.0746958) q[1];
sx q[1];
rz(-1.4935378) q[1];
sx q[1];
rz(2.1287207) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.430535) q[0];
sx q[0];
rz(-2.137888) q[0];
sx q[0];
rz(-0.60171577) q[0];
rz(-pi) q[1];
rz(2.9628721) q[2];
sx q[2];
rz(-1.0105437) q[2];
sx q[2];
rz(-0.58748875) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0433181) q[1];
sx q[1];
rz(-1.2407082) q[1];
sx q[1];
rz(1.2374452) q[1];
rz(-pi) q[2];
rz(0.49533923) q[3];
sx q[3];
rz(-2.3916349) q[3];
sx q[3];
rz(-1.0487674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0248727) q[2];
sx q[2];
rz(-2.1634384) q[2];
sx q[2];
rz(-2.8187974) q[2];
rz(2.5358477) q[3];
sx q[3];
rz(-0.79234684) q[3];
sx q[3];
rz(-0.34887031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(3.087273) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(-3.1291381) q[0];
rz(-0.74673486) q[1];
sx q[1];
rz(-2.2181999) q[1];
sx q[1];
rz(-2.8616203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075172193) q[0];
sx q[0];
rz(-0.1233347) q[0];
sx q[0];
rz(1.146011) q[0];
rz(-pi) q[1];
rz(2.0048281) q[2];
sx q[2];
rz(-1.869259) q[2];
sx q[2];
rz(-3.0500183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97023958) q[1];
sx q[1];
rz(-0.84988028) q[1];
sx q[1];
rz(-1.9848787) q[1];
rz(-pi) q[2];
rz(-1.2081228) q[3];
sx q[3];
rz(-1.8020523) q[3];
sx q[3];
rz(-1.4303007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3296457) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(-2.9296056) q[2];
rz(-0.82344615) q[3];
sx q[3];
rz(-1.4444838) q[3];
sx q[3];
rz(-2.8823891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14281808) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(-0.69277358) q[0];
rz(2.5686) q[1];
sx q[1];
rz(-1.7968105) q[1];
sx q[1];
rz(2.7105892) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26769022) q[0];
sx q[0];
rz(-2.5866716) q[0];
sx q[0];
rz(3.0303427) q[0];
rz(2.0153322) q[2];
sx q[2];
rz(-2.1241786) q[2];
sx q[2];
rz(-1.044342) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0242053) q[1];
sx q[1];
rz(-0.75953249) q[1];
sx q[1];
rz(2.7717436) q[1];
rz(-2.836976) q[3];
sx q[3];
rz(-2.0999137) q[3];
sx q[3];
rz(1.9068789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7196322) q[2];
sx q[2];
rz(-0.27406359) q[2];
sx q[2];
rz(-0.56023041) q[2];
rz(2.6719921) q[3];
sx q[3];
rz(-0.40265366) q[3];
sx q[3];
rz(-2.4277021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568759) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(-2.3003385) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(-2.2711783) q[2];
sx q[2];
rz(-0.46903789) q[2];
sx q[2];
rz(2.9901193) q[2];
rz(3.1145949) q[3];
sx q[3];
rz(-0.91201966) q[3];
sx q[3];
rz(-2.1569679) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
