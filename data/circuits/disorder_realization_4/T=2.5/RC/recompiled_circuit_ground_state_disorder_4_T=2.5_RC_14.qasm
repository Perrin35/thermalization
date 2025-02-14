OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5623915) q[0];
sx q[0];
rz(-2.3810823) q[0];
sx q[0];
rz(2.1265246) q[0];
rz(1.9782344) q[1];
sx q[1];
rz(-0.42127633) q[1];
sx q[1];
rz(1.2383229) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8684611) q[0];
sx q[0];
rz(-2.1885311) q[0];
sx q[0];
rz(0.6427838) q[0];
rz(1.1467491) q[2];
sx q[2];
rz(-1.9081665) q[2];
sx q[2];
rz(0.81126311) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1972292) q[1];
sx q[1];
rz(-2.5355593) q[1];
sx q[1];
rz(0.89971772) q[1];
rz(-1.7389033) q[3];
sx q[3];
rz(-1.013275) q[3];
sx q[3];
rz(1.3489189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6282661) q[2];
sx q[2];
rz(-1.790975) q[2];
sx q[2];
rz(2.4386621) q[2];
rz(-2.2859196) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(-1.2969016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-1.4502451) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(2.3968089) q[0];
rz(1.567747) q[1];
sx q[1];
rz(-2.4796922) q[1];
sx q[1];
rz(-2.1994798) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64025338) q[0];
sx q[0];
rz(-0.82385175) q[0];
sx q[0];
rz(-1.2938611) q[0];
x q[1];
rz(-1.6555384) q[2];
sx q[2];
rz(-2.3872445) q[2];
sx q[2];
rz(-2.1909383) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86650634) q[1];
sx q[1];
rz(-0.77436354) q[1];
sx q[1];
rz(-1.3759717) q[1];
x q[2];
rz(2.5553988) q[3];
sx q[3];
rz(-1.5206511) q[3];
sx q[3];
rz(1.8184838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.88910237) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-2.0460184) q[2];
rz(1.6628294) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(-1.3862632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(3.0304994) q[0];
sx q[0];
rz(-1.7166623) q[0];
sx q[0];
rz(-1.7362562) q[0];
rz(-1.7449215) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(-2.2415846) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42381313) q[0];
sx q[0];
rz(-2.2909145) q[0];
sx q[0];
rz(2.6507342) q[0];
rz(0.58061231) q[2];
sx q[2];
rz(-2.0502649) q[2];
sx q[2];
rz(-1.978385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.450728) q[1];
sx q[1];
rz(-0.83800661) q[1];
sx q[1];
rz(1.3370675) q[1];
x q[2];
rz(-0.95696394) q[3];
sx q[3];
rz(-1.004289) q[3];
sx q[3];
rz(1.5530128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6806543) q[2];
sx q[2];
rz(-2.9260981) q[2];
sx q[2];
rz(-2.1920965) q[2];
rz(-0.62120581) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(-1.1533823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7200274) q[0];
sx q[0];
rz(-1.255144) q[0];
sx q[0];
rz(-0.048728745) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-0.33165926) q[1];
sx q[1];
rz(-2.8299832) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6831419) q[0];
sx q[0];
rz(-1.299017) q[0];
sx q[0];
rz(-2.0608451) q[0];
rz(-pi) q[1];
rz(2.0987058) q[2];
sx q[2];
rz(-1.2706869) q[2];
sx q[2];
rz(-1.2712511) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8383465) q[1];
sx q[1];
rz(-0.7184808) q[1];
sx q[1];
rz(1.9796728) q[1];
x q[2];
rz(0.99589234) q[3];
sx q[3];
rz(-2.6290335) q[3];
sx q[3];
rz(-0.14581524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9638046) q[2];
sx q[2];
rz(-2.9298941) q[2];
sx q[2];
rz(1.4908028) q[2];
rz(2.7866411) q[3];
sx q[3];
rz(-1.4070516) q[3];
sx q[3];
rz(-1.2995592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0896924) q[0];
sx q[0];
rz(-0.35744748) q[0];
sx q[0];
rz(1.8748913) q[0];
rz(0.1162687) q[1];
sx q[1];
rz(-2.1513042) q[1];
sx q[1];
rz(1.6273392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86858803) q[0];
sx q[0];
rz(-0.75010651) q[0];
sx q[0];
rz(-2.0387893) q[0];
x q[1];
rz(1.1852988) q[2];
sx q[2];
rz(-0.23699871) q[2];
sx q[2];
rz(-1.0233819) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5942638) q[1];
sx q[1];
rz(-2.1143066) q[1];
sx q[1];
rz(-0.059652358) q[1];
x q[2];
rz(1.6125525) q[3];
sx q[3];
rz(-2.4126518) q[3];
sx q[3];
rz(1.6870013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2221471) q[2];
sx q[2];
rz(-2.6667892) q[2];
sx q[2];
rz(-1.1754645) q[2];
rz(0.90421024) q[3];
sx q[3];
rz(-1.892482) q[3];
sx q[3];
rz(2.1632532) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670052) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(2.7401155) q[0];
rz(-1.1270771) q[1];
sx q[1];
rz(-1.3104442) q[1];
sx q[1];
rz(0.81261596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75343392) q[0];
sx q[0];
rz(-1.5749314) q[0];
sx q[0];
rz(0.93664108) q[0];
x q[1];
rz(-0.27490669) q[2];
sx q[2];
rz(-1.7961575) q[2];
sx q[2];
rz(1.7079086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1620096) q[1];
sx q[1];
rz(-2.1623731) q[1];
sx q[1];
rz(1.3237557) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98315696) q[3];
sx q[3];
rz(-2.2256089) q[3];
sx q[3];
rz(-2.2838796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9408985) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(1.5852488) q[2];
rz(-0.19041348) q[3];
sx q[3];
rz(-2.3868581) q[3];
sx q[3];
rz(-1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82454005) q[0];
sx q[0];
rz(-2.7523478) q[0];
sx q[0];
rz(-0.12271605) q[0];
rz(1.9484733) q[1];
sx q[1];
rz(-2.2331608) q[1];
sx q[1];
rz(0.76748031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7372663) q[0];
sx q[0];
rz(-2.4929059) q[0];
sx q[0];
rz(-1.4859597) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38184719) q[2];
sx q[2];
rz(-1.4918461) q[2];
sx q[2];
rz(0.71982924) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0280992) q[1];
sx q[1];
rz(-0.65526795) q[1];
sx q[1];
rz(1.4987618) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010482739) q[3];
sx q[3];
rz(-0.40626486) q[3];
sx q[3];
rz(-2.3113113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7439338) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(-1.5896612) q[2];
rz(-1.6520366) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(-1.1154729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81812304) q[0];
sx q[0];
rz(-0.36190811) q[0];
sx q[0];
rz(2.7662011) q[0];
rz(2.2581532) q[1];
sx q[1];
rz(-1.6049623) q[1];
sx q[1];
rz(2.4129131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4010854) q[0];
sx q[0];
rz(-0.87941636) q[0];
sx q[0];
rz(-2.2702433) q[0];
x q[1];
rz(0.61644543) q[2];
sx q[2];
rz(-0.28646506) q[2];
sx q[2];
rz(1.1760515) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77433005) q[1];
sx q[1];
rz(-0.75241295) q[1];
sx q[1];
rz(2.1062327) q[1];
rz(-pi) q[2];
rz(0.75499423) q[3];
sx q[3];
rz(-1.0858616) q[3];
sx q[3];
rz(-0.86673966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6931307) q[2];
sx q[2];
rz(-1.5631661) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(-1.735431) q[3];
sx q[3];
rz(-1.2179255) q[3];
sx q[3];
rz(-1.5555443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2324227) q[0];
sx q[0];
rz(-0.94038832) q[0];
sx q[0];
rz(0.7255834) q[0];
rz(-1.3369417) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(-2.5673089) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5919898) q[0];
sx q[0];
rz(-1.6832325) q[0];
sx q[0];
rz(-1.7312584) q[0];
rz(-pi) q[1];
rz(-2.1404187) q[2];
sx q[2];
rz(-1.2204079) q[2];
sx q[2];
rz(0.97942615) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0813539) q[1];
sx q[1];
rz(-1.684184) q[1];
sx q[1];
rz(-2.3008623) q[1];
rz(-pi) q[2];
rz(0.78403715) q[3];
sx q[3];
rz(-1.136354) q[3];
sx q[3];
rz(-1.2617574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1672704) q[2];
sx q[2];
rz(-1.7631301) q[2];
sx q[2];
rz(-0.66217011) q[2];
rz(2.3769489) q[3];
sx q[3];
rz(-1.5352826) q[3];
sx q[3];
rz(1.3242599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26226703) q[0];
sx q[0];
rz(-0.58795324) q[0];
sx q[0];
rz(1.7437438) q[0];
rz(2.0715879) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(1.8096583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4313244) q[0];
sx q[0];
rz(-1.5402176) q[0];
sx q[0];
rz(1.4200743) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0831624) q[2];
sx q[2];
rz(-2.67423) q[2];
sx q[2];
rz(2.1455554) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3146914) q[1];
sx q[1];
rz(-1.3413635) q[1];
sx q[1];
rz(-1.8073842) q[1];
rz(-pi) q[2];
rz(0.87498192) q[3];
sx q[3];
rz(-1.1036281) q[3];
sx q[3];
rz(-0.22091978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8269044) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(2.9617214) q[2];
rz(-1.3966857) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(-2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3708645) q[0];
sx q[0];
rz(-2.6612119) q[0];
sx q[0];
rz(-2.2330855) q[0];
rz(0.52275672) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(0.54900563) q[2];
sx q[2];
rz(-0.71239757) q[2];
sx q[2];
rz(-3.0511643) q[2];
rz(2.497429) q[3];
sx q[3];
rz(-2.5732891) q[3];
sx q[3];
rz(2.7960232) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
