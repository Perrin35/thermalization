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
rz(-0.21021357) q[0];
sx q[0];
rz(-2.7855594) q[0];
sx q[0];
rz(1.5176679) q[0];
rz(1.8310945) q[1];
sx q[1];
rz(-2.2902391) q[1];
sx q[1];
rz(-0.89827615) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0097281) q[0];
sx q[0];
rz(-2.4180331) q[0];
sx q[0];
rz(2.1267164) q[0];
rz(-pi) q[1];
rz(0.29062985) q[2];
sx q[2];
rz(-2.0465474) q[2];
sx q[2];
rz(1.0889183) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3979857) q[1];
sx q[1];
rz(-1.4915257) q[1];
sx q[1];
rz(-2.958667) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4066694) q[3];
sx q[3];
rz(-1.5792819) q[3];
sx q[3];
rz(0.9134021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5400759) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(-0.69851056) q[2];
rz(1.6554333) q[3];
sx q[3];
rz(-0.58659068) q[3];
sx q[3];
rz(1.9894039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0965213) q[0];
sx q[0];
rz(-2.4424545) q[0];
sx q[0];
rz(0.29749468) q[0];
rz(2.7104764) q[1];
sx q[1];
rz(-2.2747206) q[1];
sx q[1];
rz(1.8284304) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0348059) q[0];
sx q[0];
rz(-1.2277831) q[0];
sx q[0];
rz(2.4136132) q[0];
rz(3.0975748) q[2];
sx q[2];
rz(-1.9761397) q[2];
sx q[2];
rz(2.8716033) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0370223) q[1];
sx q[1];
rz(-0.79739075) q[1];
sx q[1];
rz(-3.131011) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8937469) q[3];
sx q[3];
rz(-1.6936667) q[3];
sx q[3];
rz(0.34832277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3678652) q[2];
sx q[2];
rz(-0.44168681) q[2];
sx q[2];
rz(0.76652491) q[2];
rz(-0.68486989) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(0.49055704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3016475) q[0];
sx q[0];
rz(-1.3293581) q[0];
sx q[0];
rz(-2.8369821) q[0];
rz(2.1412762) q[1];
sx q[1];
rz(-1.1023003) q[1];
sx q[1];
rz(-2.2432378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2294126) q[0];
sx q[0];
rz(-2.1026975) q[0];
sx q[0];
rz(-0.096166111) q[0];
x q[1];
rz(-2.2730973) q[2];
sx q[2];
rz(-1.8284594) q[2];
sx q[2];
rz(-2.9106377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5542986) q[1];
sx q[1];
rz(-2.0600494) q[1];
sx q[1];
rz(-0.6502519) q[1];
x q[2];
rz(2.0933843) q[3];
sx q[3];
rz(-1.4915803) q[3];
sx q[3];
rz(-1.0635215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0072713) q[2];
sx q[2];
rz(-1.7223225) q[2];
sx q[2];
rz(-0.32709861) q[2];
rz(-1.6359811) q[3];
sx q[3];
rz(-1.4666731) q[3];
sx q[3];
rz(-1.8065642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1148249) q[0];
sx q[0];
rz(-2.6752495) q[0];
sx q[0];
rz(-2.603671) q[0];
rz(0.7750569) q[1];
sx q[1];
rz(-1.331295) q[1];
sx q[1];
rz(0.34781003) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2111172) q[0];
sx q[0];
rz(-1.0760835) q[0];
sx q[0];
rz(1.0388264) q[0];
x q[1];
rz(0.10536449) q[2];
sx q[2];
rz(-1.4658827) q[2];
sx q[2];
rz(-0.22164574) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.9197993) q[1];
sx q[1];
rz(-1.9927653) q[1];
sx q[1];
rz(-0.74733644) q[1];
rz(-pi) q[2];
rz(-0.17873879) q[3];
sx q[3];
rz(-1.1832269) q[3];
sx q[3];
rz(-2.3646858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20370087) q[2];
sx q[2];
rz(-0.72526473) q[2];
sx q[2];
rz(-2.4791278) q[2];
rz(-1.0558111) q[3];
sx q[3];
rz(-0.31848389) q[3];
sx q[3];
rz(1.8752347) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98589677) q[0];
sx q[0];
rz(-2.6031384) q[0];
sx q[0];
rz(-2.1180617) q[0];
rz(1.0515949) q[1];
sx q[1];
rz(-1.4902427) q[1];
sx q[1];
rz(-3.1255299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2002821) q[0];
sx q[0];
rz(-0.56277983) q[0];
sx q[0];
rz(-1.826123) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4532796) q[2];
sx q[2];
rz(-2.0891902) q[2];
sx q[2];
rz(-0.61358085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98707572) q[1];
sx q[1];
rz(-0.57310402) q[1];
sx q[1];
rz(0.41651233) q[1];
rz(1.315867) q[3];
sx q[3];
rz(-0.28009096) q[3];
sx q[3];
rz(-3.090429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2992799) q[2];
sx q[2];
rz(-2.2613342) q[2];
sx q[2];
rz(2.9006145) q[2];
rz(-1.109451) q[3];
sx q[3];
rz(-1.8532608) q[3];
sx q[3];
rz(2.747005) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0146765) q[0];
sx q[0];
rz(-0.83634818) q[0];
sx q[0];
rz(2.8447004) q[0];
rz(1.9706005) q[1];
sx q[1];
rz(-1.820727) q[1];
sx q[1];
rz(2.6062633) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9292262) q[0];
sx q[0];
rz(-1.3844212) q[0];
sx q[0];
rz(1.7586238) q[0];
x q[1];
rz(-2.8487327) q[2];
sx q[2];
rz(-0.83838481) q[2];
sx q[2];
rz(2.5590961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5018834) q[1];
sx q[1];
rz(-1.090126) q[1];
sx q[1];
rz(-0.23940803) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13958474) q[3];
sx q[3];
rz(-1.727316) q[3];
sx q[3];
rz(2.8536882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8239596) q[2];
sx q[2];
rz(-2.3667658) q[2];
sx q[2];
rz(-2.5114457) q[2];
rz(-1.0050425) q[3];
sx q[3];
rz(-0.21868394) q[3];
sx q[3];
rz(-2.4439243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.180535) q[0];
sx q[0];
rz(-0.068004161) q[0];
sx q[0];
rz(1.2766174) q[0];
rz(-1.0446154) q[1];
sx q[1];
rz(-1.4199384) q[1];
sx q[1];
rz(-2.9081664) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2058207) q[0];
sx q[0];
rz(-2.290831) q[0];
sx q[0];
rz(3.0854129) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4996148) q[2];
sx q[2];
rz(-1.3228088) q[2];
sx q[2];
rz(-0.61680142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9809493) q[1];
sx q[1];
rz(-2.3180974) q[1];
sx q[1];
rz(1.6232729) q[1];
rz(0.57817187) q[3];
sx q[3];
rz(-1.5147746) q[3];
sx q[3];
rz(3.0471281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.476568) q[2];
sx q[2];
rz(-1.2580405) q[2];
sx q[2];
rz(2.8705719) q[2];
rz(0.34936163) q[3];
sx q[3];
rz(-0.91785279) q[3];
sx q[3];
rz(-2.5640986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9409598) q[0];
sx q[0];
rz(-0.9372434) q[0];
sx q[0];
rz(1.8120793) q[0];
rz(-2.8726574) q[1];
sx q[1];
rz(-2.4766998) q[1];
sx q[1];
rz(1.7609133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5912) q[0];
sx q[0];
rz(-1.7730646) q[0];
sx q[0];
rz(1.3161385) q[0];
rz(-pi) q[1];
rz(2.5455238) q[2];
sx q[2];
rz(-1.5847932) q[2];
sx q[2];
rz(2.4618741) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8963651) q[1];
sx q[1];
rz(-0.62313634) q[1];
sx q[1];
rz(-1.289403) q[1];
x q[2];
rz(0.56048067) q[3];
sx q[3];
rz(-1.6874325) q[3];
sx q[3];
rz(2.9612535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5844476) q[2];
sx q[2];
rz(-2.2717805) q[2];
sx q[2];
rz(2.4073041) q[2];
rz(-1.0555438) q[3];
sx q[3];
rz(-1.5922092) q[3];
sx q[3];
rz(-2.1671104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771117) q[0];
sx q[0];
rz(-2.8900914) q[0];
sx q[0];
rz(0.86790458) q[0];
rz(-1.8382629) q[1];
sx q[1];
rz(-1.5497327) q[1];
sx q[1];
rz(-2.910639) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1385013) q[0];
sx q[0];
rz(-1.0338817) q[0];
sx q[0];
rz(2.0559539) q[0];
x q[1];
rz(-0.43416331) q[2];
sx q[2];
rz(-1.6953529) q[2];
sx q[2];
rz(-1.9366154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8924397) q[1];
sx q[1];
rz(-2.0241535) q[1];
sx q[1];
rz(-1.3282177) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1400217) q[3];
sx q[3];
rz(-2.4695307) q[3];
sx q[3];
rz(1.6035969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2086198) q[2];
sx q[2];
rz(-1.4459556) q[2];
sx q[2];
rz(-2.7033973) q[2];
rz(-2.4402601) q[3];
sx q[3];
rz(-2.0739136) q[3];
sx q[3];
rz(-2.7856538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0375131) q[0];
sx q[0];
rz(-0.47261819) q[0];
sx q[0];
rz(-2.0613101) q[0];
rz(2.089962) q[1];
sx q[1];
rz(-2.0104505) q[1];
sx q[1];
rz(1.0952605) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0290463) q[0];
sx q[0];
rz(-1.8575107) q[0];
sx q[0];
rz(-2.6081144) q[0];
rz(-pi) q[1];
rz(3.0313086) q[2];
sx q[2];
rz(-0.37026893) q[2];
sx q[2];
rz(-0.25843378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.238823) q[1];
sx q[1];
rz(-2.2939957) q[1];
sx q[1];
rz(-0.1478424) q[1];
rz(1.0686103) q[3];
sx q[3];
rz(-1.0334224) q[3];
sx q[3];
rz(0.93056941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1209391) q[2];
sx q[2];
rz(-1.8755269) q[2];
sx q[2];
rz(0.90547639) q[2];
rz(-0.71546537) q[3];
sx q[3];
rz(-0.72532907) q[3];
sx q[3];
rz(-0.65413094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41659551) q[0];
sx q[0];
rz(-1.409197) q[0];
sx q[0];
rz(-2.5247164) q[0];
rz(-1.1299409) q[1];
sx q[1];
rz(-1.8604953) q[1];
sx q[1];
rz(-1.4166191) q[1];
rz(2.2466725) q[2];
sx q[2];
rz(-0.85325675) q[2];
sx q[2];
rz(0.85430145) q[2];
rz(1.7814287) q[3];
sx q[3];
rz(-1.7192626) q[3];
sx q[3];
rz(0.33178793) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
