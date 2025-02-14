OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4717785) q[0];
sx q[0];
rz(-1.9999003) q[0];
sx q[0];
rz(0.25622955) q[0];
rz(3.5186634) q[1];
sx q[1];
rz(4.9405603) q[1];
sx q[1];
rz(8.3648051) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621583) q[0];
sx q[0];
rz(-0.16955626) q[0];
sx q[0];
rz(2.3441699) q[0];
rz(-2.4490812) q[2];
sx q[2];
rz(-1.5438809) q[2];
sx q[2];
rz(-1.4262425) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.022885325) q[1];
sx q[1];
rz(-1.361025) q[1];
sx q[1];
rz(-0.24521141) q[1];
x q[2];
rz(2.9971231) q[3];
sx q[3];
rz(-1.472755) q[3];
sx q[3];
rz(-1.4812183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18067351) q[2];
sx q[2];
rz(-0.91327614) q[2];
sx q[2];
rz(-1.712435) q[2];
rz(2.5299431) q[3];
sx q[3];
rz(-1.6983756) q[3];
sx q[3];
rz(0.071368607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3058158) q[0];
sx q[0];
rz(-0.80128765) q[0];
sx q[0];
rz(-0.74547705) q[0];
rz(-1.2019134) q[1];
sx q[1];
rz(-1.3246091) q[1];
sx q[1];
rz(2.1048996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03464493) q[0];
sx q[0];
rz(-1.4781524) q[0];
sx q[0];
rz(-1.205446) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5971423) q[2];
sx q[2];
rz(-2.3374632) q[2];
sx q[2];
rz(-2.1955817) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.17192852) q[1];
sx q[1];
rz(-1.2873889) q[1];
sx q[1];
rz(0.71782459) q[1];
rz(-pi) q[2];
rz(-0.3221945) q[3];
sx q[3];
rz(-0.60880843) q[3];
sx q[3];
rz(-2.6170066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1408954) q[2];
sx q[2];
rz(-1.2709728) q[2];
sx q[2];
rz(-0.97266436) q[2];
rz(2.2843212) q[3];
sx q[3];
rz(-2.075115) q[3];
sx q[3];
rz(1.8734141) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.641441) q[0];
sx q[0];
rz(-0.86857906) q[0];
sx q[0];
rz(0.80297536) q[0];
rz(-2.4909486) q[1];
sx q[1];
rz(-0.79901189) q[1];
sx q[1];
rz(2.9237936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9150216) q[0];
sx q[0];
rz(-2.3816097) q[0];
sx q[0];
rz(1.0081916) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8286092) q[2];
sx q[2];
rz(-0.65201983) q[2];
sx q[2];
rz(-2.5283732) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4068702) q[1];
sx q[1];
rz(-2.0592505) q[1];
sx q[1];
rz(3.0367756) q[1];
x q[2];
rz(-0.33588107) q[3];
sx q[3];
rz(-1.7917197) q[3];
sx q[3];
rz(0.22798746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.190072) q[2];
sx q[2];
rz(-2.0900487) q[2];
sx q[2];
rz(-1.5104843) q[2];
rz(-1.2269646) q[3];
sx q[3];
rz(-1.0534143) q[3];
sx q[3];
rz(2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25201061) q[0];
sx q[0];
rz(-2.4376106) q[0];
sx q[0];
rz(-2.4776283) q[0];
rz(0.12531677) q[1];
sx q[1];
rz(-1.434451) q[1];
sx q[1];
rz(-2.1024316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0200203) q[0];
sx q[0];
rz(-1.5680321) q[0];
sx q[0];
rz(0.0098258709) q[0];
x q[1];
rz(2.8739615) q[2];
sx q[2];
rz(-1.1881141) q[2];
sx q[2];
rz(-3.1185993) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.127847) q[1];
sx q[1];
rz(-1.7798073) q[1];
sx q[1];
rz(-2.0621939) q[1];
x q[2];
rz(-2.6986946) q[3];
sx q[3];
rz(-1.6567461) q[3];
sx q[3];
rz(-1.3493054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.30109721) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(3.0636129) q[2];
rz(1.1178499) q[3];
sx q[3];
rz(-2.100914) q[3];
sx q[3];
rz(1.0361205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(0.41314405) q[0];
rz(-1.9059937) q[1];
sx q[1];
rz(-1.3256336) q[1];
sx q[1];
rz(1.3135757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298656) q[0];
sx q[0];
rz(-2.4851375) q[0];
sx q[0];
rz(1.2675499) q[0];
rz(-2.0356455) q[2];
sx q[2];
rz(-1.8951956) q[2];
sx q[2];
rz(1.731763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2755679) q[1];
sx q[1];
rz(-1.1236262) q[1];
sx q[1];
rz(-2.7204896) q[1];
x q[2];
rz(-2.026346) q[3];
sx q[3];
rz(-0.83720647) q[3];
sx q[3];
rz(2.7546492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95973394) q[2];
sx q[2];
rz(-1.2343957) q[2];
sx q[2];
rz(-2.6971297) q[2];
rz(-1.9410761) q[3];
sx q[3];
rz(-1.2512755) q[3];
sx q[3];
rz(-3.0357231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34894094) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(-3.0643903) q[0];
rz(2.1791747) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(-0.033800689) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0173967) q[0];
sx q[0];
rz(-0.58308812) q[0];
sx q[0];
rz(0.65753196) q[0];
rz(3.097418) q[2];
sx q[2];
rz(-0.54810537) q[2];
sx q[2];
rz(-1.7063315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56735605) q[1];
sx q[1];
rz(-2.1891382) q[1];
sx q[1];
rz(0.38783698) q[1];
rz(-pi) q[2];
rz(-1.46557) q[3];
sx q[3];
rz(-0.79685539) q[3];
sx q[3];
rz(-1.5060671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17322156) q[2];
sx q[2];
rz(-1.6049478) q[2];
sx q[2];
rz(2.8823493) q[2];
rz(0.94414583) q[3];
sx q[3];
rz(-0.21374948) q[3];
sx q[3];
rz(1.8102144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060870085) q[0];
sx q[0];
rz(-0.20651564) q[0];
sx q[0];
rz(-0.61332214) q[0];
rz(-2.2382286) q[1];
sx q[1];
rz(-2.3009243) q[1];
sx q[1];
rz(-0.14796743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4658303) q[0];
sx q[0];
rz(-2.0480541) q[0];
sx q[0];
rz(-2.3289519) q[0];
rz(-pi) q[1];
rz(0.067909165) q[2];
sx q[2];
rz(-0.43988827) q[2];
sx q[2];
rz(1.2665757) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92098713) q[1];
sx q[1];
rz(-1.9927532) q[1];
sx q[1];
rz(-0.76442952) q[1];
x q[2];
rz(1.5336125) q[3];
sx q[3];
rz(-2.5761009) q[3];
sx q[3];
rz(2.9786106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.938574) q[2];
sx q[2];
rz(-1.729894) q[2];
sx q[2];
rz(2.1417248) q[2];
rz(1.3153007) q[3];
sx q[3];
rz(-0.2838997) q[3];
sx q[3];
rz(0.6704754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.1891747) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(2.1536105) q[0];
rz(1.8723764) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(0.16407897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7456147) q[0];
sx q[0];
rz(-1.656129) q[0];
sx q[0];
rz(3.0239168) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0751749) q[2];
sx q[2];
rz(-1.7640503) q[2];
sx q[2];
rz(2.9444864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33262256) q[1];
sx q[1];
rz(-1.5738942) q[1];
sx q[1];
rz(-1.078152) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2211391) q[3];
sx q[3];
rz(-1.4152539) q[3];
sx q[3];
rz(2.807164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24171955) q[2];
sx q[2];
rz(-2.6417929) q[2];
sx q[2];
rz(2.0797753) q[2];
rz(-0.1344943) q[3];
sx q[3];
rz(-0.89956346) q[3];
sx q[3];
rz(-0.053248052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6580842) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(-0.2051951) q[0];
rz(2.9070053) q[1];
sx q[1];
rz(-1.4261475) q[1];
sx q[1];
rz(2.9451784) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0935555) q[0];
sx q[0];
rz(-1.8797848) q[0];
sx q[0];
rz(0.97889203) q[0];
rz(0.5242879) q[2];
sx q[2];
rz(-1.4349738) q[2];
sx q[2];
rz(0.4834396) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82746303) q[1];
sx q[1];
rz(-1.8937096) q[1];
sx q[1];
rz(0.18675487) q[1];
x q[2];
rz(-2.6077723) q[3];
sx q[3];
rz(-2.2861423) q[3];
sx q[3];
rz(1.2792339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3966169) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(1.9632001) q[2];
rz(0.34934238) q[3];
sx q[3];
rz(-2.0125466) q[3];
sx q[3];
rz(-1.0497302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1402533) q[0];
sx q[0];
rz(-1.1118735) q[0];
sx q[0];
rz(2.4358791) q[0];
rz(0.69752518) q[1];
sx q[1];
rz(-1.3730647) q[1];
sx q[1];
rz(-0.90550214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5030497) q[0];
sx q[0];
rz(-2.7075504) q[0];
sx q[0];
rz(0.040157138) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72737965) q[2];
sx q[2];
rz(-0.47107163) q[2];
sx q[2];
rz(-1.7076275) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9230629) q[1];
sx q[1];
rz(-2.8273812) q[1];
sx q[1];
rz(-0.15337069) q[1];
rz(-pi) q[2];
rz(-1.2071183) q[3];
sx q[3];
rz(-1.26059) q[3];
sx q[3];
rz(0.23387533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2424348) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(-0.71851292) q[2];
rz(2.4588623) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(-0.1040641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.318442) q[0];
sx q[0];
rz(-0.57118509) q[0];
sx q[0];
rz(-0.1027064) q[0];
rz(1.5009343) q[1];
sx q[1];
rz(-1.2702912) q[1];
sx q[1];
rz(2.695695) q[1];
rz(-1.6886383) q[2];
sx q[2];
rz(-1.4284882) q[2];
sx q[2];
rz(1.535648) q[2];
rz(2.4950776) q[3];
sx q[3];
rz(-2.2263891) q[3];
sx q[3];
rz(1.9937766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
