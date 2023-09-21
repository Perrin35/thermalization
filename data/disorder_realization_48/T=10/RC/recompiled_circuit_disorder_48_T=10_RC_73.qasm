OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(3.9324023) q[0];
sx q[0];
rz(12.232236) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96555644) q[0];
sx q[0];
rz(-0.82057014) q[0];
sx q[0];
rz(-1.6198938) q[0];
rz(1.5300418) q[2];
sx q[2];
rz(-1.1062804) q[2];
sx q[2];
rz(-1.0711311) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3764066) q[1];
sx q[1];
rz(-1.1131439) q[1];
sx q[1];
rz(-2.5916369) q[1];
rz(-1.3487885) q[3];
sx q[3];
rz(-1.3862002) q[3];
sx q[3];
rz(0.27664646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.618764) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(-0.5775601) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(-2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(-1.739025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4336006) q[0];
sx q[0];
rz(-0.029768243) q[0];
sx q[0];
rz(3.004651) q[0];
rz(-pi) q[1];
rz(0.34105532) q[2];
sx q[2];
rz(-0.58652069) q[2];
sx q[2];
rz(-1.7797433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.95784159) q[1];
sx q[1];
rz(-1.9722003) q[1];
sx q[1];
rz(-2.6028231) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6321294) q[3];
sx q[3];
rz(-1.4940133) q[3];
sx q[3];
rz(2.0413105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7188321) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3690255) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(2.7405222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7558407) q[0];
sx q[0];
rz(-0.34163654) q[0];
sx q[0];
rz(-1.2582448) q[0];
x q[1];
rz(0.24984078) q[2];
sx q[2];
rz(-0.65953883) q[2];
sx q[2];
rz(-1.8674873) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9123147) q[1];
sx q[1];
rz(-2.4114128) q[1];
sx q[1];
rz(-3.1320523) q[1];
x q[2];
rz(-1.1776351) q[3];
sx q[3];
rz(-1.7522246) q[3];
sx q[3];
rz(1.7593311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(-2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.697631) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(-2.8881883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4454173) q[0];
sx q[0];
rz(-1.6735958) q[0];
sx q[0];
rz(-1.0284543) q[0];
x q[1];
rz(0.01732145) q[2];
sx q[2];
rz(-2.0565196) q[2];
sx q[2];
rz(-2.3522365) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.464444) q[1];
sx q[1];
rz(-2.1975937) q[1];
sx q[1];
rz(1.6897175) q[1];
rz(2.4921791) q[3];
sx q[3];
rz(-1.489349) q[3];
sx q[3];
rz(-2.7922975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0306586) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(0.91286719) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(-3.0786247) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(0.45809349) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23819085) q[0];
sx q[0];
rz(-1.5731249) q[0];
sx q[0];
rz(2.6989614) q[0];
x q[1];
rz(-0.45256726) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(2.3134311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8372247) q[1];
sx q[1];
rz(-0.75290426) q[1];
sx q[1];
rz(1.3992845) q[1];
x q[2];
rz(0.094893806) q[3];
sx q[3];
rz(-1.0720383) q[3];
sx q[3];
rz(2.9961078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(-3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620419) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(-1.8744291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8454682) q[0];
sx q[0];
rz(-2.2757029) q[0];
sx q[0];
rz(2.092093) q[0];
x q[1];
rz(2.3976372) q[2];
sx q[2];
rz(-1.4641054) q[2];
sx q[2];
rz(-2.8731186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8900745) q[1];
sx q[1];
rz(-1.8971844) q[1];
sx q[1];
rz(-0.20784394) q[1];
rz(-pi) q[2];
x q[2];
rz(1.091286) q[3];
sx q[3];
rz(-1.7877842) q[3];
sx q[3];
rz(-3.019141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(2.069058) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(-1.1118719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8060914) q[0];
sx q[0];
rz(-1.6643235) q[0];
sx q[0];
rz(1.085824) q[0];
rz(-pi) q[1];
rz(-2.6107437) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(-1.2209148) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.062473) q[1];
sx q[1];
rz(-1.2812496) q[1];
sx q[1];
rz(1.6505961) q[1];
x q[2];
rz(-2.7518919) q[3];
sx q[3];
rz(-1.1003564) q[3];
sx q[3];
rz(2.5999992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83773461) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(1.6052823) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.4060098) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(-0.73079601) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(-2.3866167) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6801493) q[0];
sx q[0];
rz(-1.7168683) q[0];
sx q[0];
rz(0.84772528) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63379143) q[2];
sx q[2];
rz(-0.83665028) q[2];
sx q[2];
rz(1.447669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.054143993) q[1];
sx q[1];
rz(-0.99522299) q[1];
sx q[1];
rz(1.4364442) q[1];
rz(-pi) q[2];
rz(-2.2488392) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(1.0196109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-0.63684741) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.4144142) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(2.0027347) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(-3.1220904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0855904) q[0];
sx q[0];
rz(-1.7770355) q[0];
sx q[0];
rz(-0.081736728) q[0];
rz(-pi) q[1];
rz(0.15140622) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(0.99934794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86917415) q[1];
sx q[1];
rz(-2.4915016) q[1];
sx q[1];
rz(1.2692578) q[1];
x q[2];
rz(1.8832302) q[3];
sx q[3];
rz(-1.7282681) q[3];
sx q[3];
rz(1.0851932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-2.2926889) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983109) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(-1.7548521) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.9932995) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5507817) q[0];
sx q[0];
rz(-1.592144) q[0];
sx q[0];
rz(-0.14069964) q[0];
rz(-0.30015595) q[2];
sx q[2];
rz(-2.488392) q[2];
sx q[2];
rz(-2.4094827) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6069968) q[1];
sx q[1];
rz(-1.7417223) q[1];
sx q[1];
rz(2.2439438) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5185353) q[3];
sx q[3];
rz(-1.9196379) q[3];
sx q[3];
rz(0.40961743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2293573) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0317595) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-0.96314349) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(2.4776496) q[2];
sx q[2];
rz(-1.2524111) q[2];
sx q[2];
rz(-1.8128957) q[2];
rz(1.6297324) q[3];
sx q[3];
rz(-2.5806576) q[3];
sx q[3];
rz(2.9604119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];