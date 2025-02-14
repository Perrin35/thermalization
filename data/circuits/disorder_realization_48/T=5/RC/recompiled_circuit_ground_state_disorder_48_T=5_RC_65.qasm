OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(-1.1416924) q[0];
sx q[0];
rz(-0.25622955) q[0];
rz(-2.7645219) q[1];
sx q[1];
rz(-1.3426251) q[1];
sx q[1];
rz(-1.0599729) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7840756) q[0];
sx q[0];
rz(-1.6889484) q[0];
sx q[0];
rz(-1.4489001) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5358309) q[2];
sx q[2];
rz(-0.87858534) q[2];
sx q[2];
rz(-0.12223211) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.022885325) q[1];
sx q[1];
rz(-1.361025) q[1];
sx q[1];
rz(-0.24521141) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5422402) q[3];
sx q[3];
rz(-0.17440344) q[3];
sx q[3];
rz(2.6389183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.18067351) q[2];
sx q[2];
rz(-0.91327614) q[2];
sx q[2];
rz(-1.4291576) q[2];
rz(2.5299431) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(3.070224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83577689) q[0];
sx q[0];
rz(-0.80128765) q[0];
sx q[0];
rz(2.3961156) q[0];
rz(-1.2019134) q[1];
sx q[1];
rz(-1.3246091) q[1];
sx q[1];
rz(2.1048996) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1069477) q[0];
sx q[0];
rz(-1.4781524) q[0];
sx q[0];
rz(1.9361467) q[0];
x q[1];
rz(1.5971423) q[2];
sx q[2];
rz(-0.80412946) q[2];
sx q[2];
rz(-0.94601099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0889609) q[1];
sx q[1];
rz(-0.76238576) q[1];
sx q[1];
rz(2.7247468) q[1];
rz(-pi) q[2];
rz(-1.3535301) q[3];
sx q[3];
rz(-2.144162) q[3];
sx q[3];
rz(2.2305302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0006973) q[2];
sx q[2];
rz(-1.2709728) q[2];
sx q[2];
rz(2.1689283) q[2];
rz(0.85727143) q[3];
sx q[3];
rz(-2.075115) q[3];
sx q[3];
rz(-1.8734141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.641441) q[0];
sx q[0];
rz(-0.86857906) q[0];
sx q[0];
rz(2.3386173) q[0];
rz(0.650644) q[1];
sx q[1];
rz(-0.79901189) q[1];
sx q[1];
rz(2.9237936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.652199) q[0];
sx q[0];
rz(-2.1928761) q[0];
sx q[0];
rz(2.6724044) q[0];
x q[1];
rz(1.8016544) q[2];
sx q[2];
rz(-0.95536026) q[2];
sx q[2];
rz(0.99986693) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1863292) q[1];
sx q[1];
rz(-2.6429089) q[1];
sx q[1];
rz(-1.3763983) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8043376) q[3];
sx q[3];
rz(-1.8982049) q[3];
sx q[3];
rz(1.722432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.190072) q[2];
sx q[2];
rz(-2.0900487) q[2];
sx q[2];
rz(1.6311084) q[2];
rz(1.914628) q[3];
sx q[3];
rz(-1.0534143) q[3];
sx q[3];
rz(-1.0043043) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.889582) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(-2.4776283) q[0];
rz(3.0162759) q[1];
sx q[1];
rz(-1.7071416) q[1];
sx q[1];
rz(-2.1024316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44925112) q[0];
sx q[0];
rz(-1.5609705) q[0];
sx q[0];
rz(-1.5735606) q[0];
x q[1];
rz(2.8739615) q[2];
sx q[2];
rz(-1.9534785) q[2];
sx q[2];
rz(3.1185993) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.127847) q[1];
sx q[1];
rz(-1.3617853) q[1];
sx q[1];
rz(1.0793988) q[1];
x q[2];
rz(1.4757206) q[3];
sx q[3];
rz(-2.0119442) q[3];
sx q[3];
rz(-2.8794062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8404954) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(0.077979716) q[2];
rz(1.1178499) q[3];
sx q[3];
rz(-1.0406787) q[3];
sx q[3];
rz(2.1054721) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2074821) q[0];
sx q[0];
rz(-0.20620646) q[0];
sx q[0];
rz(-2.7284486) q[0];
rz(-1.9059937) q[1];
sx q[1];
rz(-1.3256336) q[1];
sx q[1];
rz(1.3135757) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020371) q[0];
sx q[0];
rz(-1.7540723) q[0];
sx q[0];
rz(0.9367783) q[0];
rz(1.1059472) q[2];
sx q[2];
rz(-1.2463971) q[2];
sx q[2];
rz(-1.731763) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6690065) q[1];
sx q[1];
rz(-0.60429497) q[1];
sx q[1];
rz(-0.86493203) q[1];
rz(-2.3543139) q[3];
sx q[3];
rz(-1.2378927) q[3];
sx q[3];
rz(1.5008139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95973394) q[2];
sx q[2];
rz(-1.9071969) q[2];
sx q[2];
rz(-0.44446298) q[2];
rz(-1.9410761) q[3];
sx q[3];
rz(-1.2512755) q[3];
sx q[3];
rz(0.1058696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926517) q[0];
sx q[0];
rz(-0.3949202) q[0];
sx q[0];
rz(3.0643903) q[0];
rz(-2.1791747) q[1];
sx q[1];
rz(-1.9541157) q[1];
sx q[1];
rz(3.107792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156436) q[0];
sx q[0];
rz(-1.2275877) q[0];
sx q[0];
rz(-0.48114602) q[0];
rz(3.097418) q[2];
sx q[2];
rz(-0.54810537) q[2];
sx q[2];
rz(1.4352611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0950468) q[1];
sx q[1];
rz(-2.425417) q[1];
sx q[1];
rz(-2.0594055) q[1];
x q[2];
rz(2.364879) q[3];
sx q[3];
rz(-1.6459822) q[3];
sx q[3];
rz(0.0089587072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9683711) q[2];
sx q[2];
rz(-1.5366448) q[2];
sx q[2];
rz(-2.8823493) q[2];
rz(-0.94414583) q[3];
sx q[3];
rz(-0.21374948) q[3];
sx q[3];
rz(1.3313782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.060870085) q[0];
sx q[0];
rz(-0.20651564) q[0];
sx q[0];
rz(2.5282705) q[0];
rz(0.90336409) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(-2.9936252) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4658303) q[0];
sx q[0];
rz(-2.0480541) q[0];
sx q[0];
rz(2.3289519) q[0];
rz(-1.602722) q[2];
sx q[2];
rz(-2.0096001) q[2];
sx q[2];
rz(1.9500458) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8954249) q[1];
sx q[1];
rz(-0.85201272) q[1];
sx q[1];
rz(-2.5661929) q[1];
rz(-pi) q[2];
rz(1.0056173) q[3];
sx q[3];
rz(-1.5508754) q[3];
sx q[3];
rz(-1.3764149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.938574) q[2];
sx q[2];
rz(-1.4116986) q[2];
sx q[2];
rz(2.1417248) q[2];
rz(-1.3153007) q[3];
sx q[3];
rz(-2.857693) q[3];
sx q[3];
rz(0.6704754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95241791) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(-2.1536105) q[0];
rz(-1.2692163) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(0.16407897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39597797) q[0];
sx q[0];
rz(-1.4854637) q[0];
sx q[0];
rz(3.0239168) q[0];
x q[1];
rz(-2.0664178) q[2];
sx q[2];
rz(-1.3775423) q[2];
sx q[2];
rz(-0.19710625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8089701) q[1];
sx q[1];
rz(-1.5676985) q[1];
sx q[1];
rz(-2.0634406) q[1];
rz(0.19453519) q[3];
sx q[3];
rz(-2.2119868) q[3];
sx q[3];
rz(-1.3536782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24171955) q[2];
sx q[2];
rz(-0.49979979) q[2];
sx q[2];
rz(1.0618173) q[2];
rz(-3.0070983) q[3];
sx q[3];
rz(-0.89956346) q[3];
sx q[3];
rz(0.053248052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6580842) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(2.9363976) q[0];
rz(0.23458734) q[1];
sx q[1];
rz(-1.4261475) q[1];
sx q[1];
rz(-2.9451784) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.043606) q[0];
sx q[0];
rz(-0.65904407) q[0];
sx q[0];
rz(1.0511257) q[0];
x q[1];
rz(-0.5242879) q[2];
sx q[2];
rz(-1.7066188) q[2];
sx q[2];
rz(0.4834396) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3383693) q[1];
sx q[1];
rz(-1.3937989) q[1];
sx q[1];
rz(-1.2425674) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78077448) q[3];
sx q[3];
rz(-1.1765683) q[3];
sx q[3];
rz(3.0633283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74497574) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(1.1783925) q[2];
rz(-0.34934238) q[3];
sx q[3];
rz(-1.1290461) q[3];
sx q[3];
rz(-1.0497302) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1402533) q[0];
sx q[0];
rz(-2.0297191) q[0];
sx q[0];
rz(-2.4358791) q[0];
rz(-0.69752518) q[1];
sx q[1];
rz(-1.768528) q[1];
sx q[1];
rz(2.2360905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10418358) q[0];
sx q[0];
rz(-1.5539123) q[0];
sx q[0];
rz(-2.707858) q[0];
rz(-2.7780813) q[2];
sx q[2];
rz(-1.8773407) q[2];
sx q[2];
rz(0.53378045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21852979) q[1];
sx q[1];
rz(-0.31421146) q[1];
sx q[1];
rz(2.988222) q[1];
x q[2];
rz(-2.3042514) q[3];
sx q[3];
rz(-2.6681134) q[3];
sx q[3];
rz(2.480643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2424348) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(-2.4230797) q[2];
rz(2.4588623) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(-0.1040641) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8231507) q[0];
sx q[0];
rz(-2.5704076) q[0];
sx q[0];
rz(3.0388863) q[0];
rz(1.6406583) q[1];
sx q[1];
rz(-1.8713015) q[1];
sx q[1];
rz(-0.4458977) q[1];
rz(-2.9983042) q[2];
sx q[2];
rz(-1.454151) q[2];
sx q[2];
rz(3.1232338) q[2];
rz(-0.90632306) q[3];
sx q[3];
rz(-2.255848) q[3];
sx q[3];
rz(-2.0391603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
