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
rz(2.9313791) q[0];
sx q[0];
rz(5.9271521) q[0];
sx q[0];
rz(7.9071101) q[0];
rz(-1.3104982) q[1];
sx q[1];
rz(-0.85135353) q[1];
sx q[1];
rz(0.89827615) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87469989) q[0];
sx q[0];
rz(-1.9277097) q[0];
sx q[0];
rz(-0.92706417) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29062985) q[2];
sx q[2];
rz(-1.0950452) q[2];
sx q[2];
rz(2.0526744) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15816244) q[1];
sx q[1];
rz(-1.753141) q[1];
sx q[1];
rz(1.4901864) q[1];
rz(-0.0086011767) q[3];
sx q[3];
rz(-1.4066753) q[3];
sx q[3];
rz(2.4827931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6015168) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(-0.69851056) q[2];
rz(-1.4861594) q[3];
sx q[3];
rz(-0.58659068) q[3];
sx q[3];
rz(1.9894039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0965213) q[0];
sx q[0];
rz(-2.4424545) q[0];
sx q[0];
rz(-0.29749468) q[0];
rz(-0.43111626) q[1];
sx q[1];
rz(-0.86687207) q[1];
sx q[1];
rz(1.3131622) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8968412) q[0];
sx q[0];
rz(-2.2479821) q[0];
sx q[0];
rz(-2.0170101) q[0];
rz(-pi) q[1];
rz(-1.976491) q[2];
sx q[2];
rz(-1.6112453) q[2];
sx q[2];
rz(1.8581529) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6152108) q[1];
sx q[1];
rz(-1.5783678) q[1];
sx q[1];
rz(-0.79736276) q[1];
x q[2];
rz(-2.8937469) q[3];
sx q[3];
rz(-1.447926) q[3];
sx q[3];
rz(-0.34832277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7737274) q[2];
sx q[2];
rz(-2.6999058) q[2];
sx q[2];
rz(-0.76652491) q[2];
rz(-2.4567228) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(2.6510356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3016475) q[0];
sx q[0];
rz(-1.8122346) q[0];
sx q[0];
rz(-0.30461052) q[0];
rz(-2.1412762) q[1];
sx q[1];
rz(-2.0392923) q[1];
sx q[1];
rz(0.89835483) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91218006) q[0];
sx q[0];
rz(-1.0388952) q[0];
sx q[0];
rz(-0.096166111) q[0];
rz(-pi) q[1];
rz(-0.33240357) q[2];
sx q[2];
rz(-0.89611182) q[2];
sx q[2];
rz(-1.5522267) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5542986) q[1];
sx q[1];
rz(-1.0815433) q[1];
sx q[1];
rz(2.4913408) q[1];
rz(-1.4130745) q[3];
sx q[3];
rz(-0.52800679) q[3];
sx q[3];
rz(0.64380336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0072713) q[2];
sx q[2];
rz(-1.4192702) q[2];
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
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026767749) q[0];
sx q[0];
rz(-2.6752495) q[0];
sx q[0];
rz(-0.53792167) q[0];
rz(-0.7750569) q[1];
sx q[1];
rz(-1.331295) q[1];
sx q[1];
rz(2.7937826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7737601) q[0];
sx q[0];
rz(-2.0335541) q[0];
sx q[0];
rz(-0.55928834) q[0];
rz(-3.0362282) q[2];
sx q[2];
rz(-1.4658827) q[2];
sx q[2];
rz(-0.22164574) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0744275) q[1];
sx q[1];
rz(-0.8377155) q[1];
sx q[1];
rz(0.58372402) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9640173) q[3];
sx q[3];
rz(-1.4054417) q[3];
sx q[3];
rz(-2.4158784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9378918) q[2];
sx q[2];
rz(-2.4163279) q[2];
sx q[2];
rz(-0.66246486) q[2];
rz(1.0558111) q[3];
sx q[3];
rz(-2.8231088) q[3];
sx q[3];
rz(-1.2663579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1556959) q[0];
sx q[0];
rz(-0.53845423) q[0];
sx q[0];
rz(-2.1180617) q[0];
rz(2.0899978) q[1];
sx q[1];
rz(-1.4902427) q[1];
sx q[1];
rz(-0.016062707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4996289) q[0];
sx q[0];
rz(-2.1132541) q[0];
sx q[0];
rz(0.15799518) q[0];
rz(-0.73170029) q[2];
sx q[2];
rz(-2.3061487) q[2];
sx q[2];
rz(-1.6421865) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2019883) q[1];
sx q[1];
rz(-1.7919722) q[1];
sx q[1];
rz(2.6084234) q[1];
rz(-pi) q[2];
rz(-1.2993126) q[3];
sx q[3];
rz(-1.6405655) q[3];
sx q[3];
rz(1.8673563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2992799) q[2];
sx q[2];
rz(-2.2613342) q[2];
sx q[2];
rz(-0.24097815) q[2];
rz(-1.109451) q[3];
sx q[3];
rz(-1.8532608) q[3];
sx q[3];
rz(-0.39458767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0146765) q[0];
sx q[0];
rz(-0.83634818) q[0];
sx q[0];
rz(2.8447004) q[0];
rz(-1.1709921) q[1];
sx q[1];
rz(-1.3208656) q[1];
sx q[1];
rz(-2.6062633) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7479582) q[0];
sx q[0];
rz(-1.7553333) q[0];
sx q[0];
rz(2.9519597) q[0];
rz(-pi) q[1];
rz(1.2601603) q[2];
sx q[2];
rz(-0.77858965) q[2];
sx q[2];
rz(2.1354577) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1254309) q[1];
sx q[1];
rz(-2.6088073) q[1];
sx q[1];
rz(1.1440117) q[1];
x q[2];
rz(1.728828) q[3];
sx q[3];
rz(-1.4329289) q[3];
sx q[3];
rz(1.3047895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.31763306) q[2];
sx q[2];
rz(-2.3667658) q[2];
sx q[2];
rz(0.63014692) q[2];
rz(-1.0050425) q[3];
sx q[3];
rz(-2.9229087) q[3];
sx q[3];
rz(-0.69766831) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.180535) q[0];
sx q[0];
rz(-0.068004161) q[0];
sx q[0];
rz(-1.2766174) q[0];
rz(2.0969773) q[1];
sx q[1];
rz(-1.4199384) q[1];
sx q[1];
rz(0.23342625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.935772) q[0];
sx q[0];
rz(-0.85076166) q[0];
sx q[0];
rz(-0.056179742) q[0];
x q[1];
rz(-0.40005513) q[2];
sx q[2];
rz(-2.4597485) q[2];
sx q[2];
rz(1.2712196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.4458369) q[1];
sx q[1];
rz(-1.6092811) q[1];
sx q[1];
rz(-0.74798788) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6376611) q[3];
sx q[3];
rz(-2.1479448) q[3];
sx q[3];
rz(1.4398097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.6650247) q[2];
sx q[2];
rz(-1.8835521) q[2];
sx q[2];
rz(-0.27102077) q[2];
rz(-2.792231) q[3];
sx q[3];
rz(-2.2237399) q[3];
sx q[3];
rz(2.5640986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9409598) q[0];
sx q[0];
rz(-2.2043493) q[0];
sx q[0];
rz(-1.8120793) q[0];
rz(2.8726574) q[1];
sx q[1];
rz(-2.4766998) q[1];
sx q[1];
rz(-1.7609133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5047245) q[0];
sx q[0];
rz(-0.32384017) q[0];
sx q[0];
rz(2.25405) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59606885) q[2];
sx q[2];
rz(-1.5847932) q[2];
sx q[2];
rz(-0.67971855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.24522752) q[1];
sx q[1];
rz(-0.62313634) q[1];
sx q[1];
rz(1.8521897) q[1];
rz(-1.4333357) q[3];
sx q[3];
rz(-1.0145742) q[3];
sx q[3];
rz(1.463365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5844476) q[2];
sx q[2];
rz(-0.86981213) q[2];
sx q[2];
rz(-0.73428854) q[2];
rz(2.0860489) q[3];
sx q[3];
rz(-1.5493834) q[3];
sx q[3];
rz(2.1671104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.664481) q[0];
sx q[0];
rz(-0.25150126) q[0];
sx q[0];
rz(2.2736881) q[0];
rz(-1.3033298) q[1];
sx q[1];
rz(-1.5918599) q[1];
sx q[1];
rz(0.23095362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.445914) q[0];
sx q[0];
rz(-1.983108) q[0];
sx q[0];
rz(-0.59230174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7079389) q[2];
sx q[2];
rz(-2.0013705) q[2];
sx q[2];
rz(0.42335934) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7119904) q[1];
sx q[1];
rz(-1.3531405) q[1];
sx q[1];
rz(2.6764286) q[1];
rz(-pi) q[2];
rz(-2.4695314) q[3];
sx q[3];
rz(-1.5717744) q[3];
sx q[3];
rz(3.1100215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2086198) q[2];
sx q[2];
rz(-1.695637) q[2];
sx q[2];
rz(-2.7033973) q[2];
rz(2.4402601) q[3];
sx q[3];
rz(-1.067679) q[3];
sx q[3];
rz(0.35593885) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1040795) q[0];
sx q[0];
rz(-0.47261819) q[0];
sx q[0];
rz(-1.0802826) q[0];
rz(-2.089962) q[1];
sx q[1];
rz(-2.0104505) q[1];
sx q[1];
rz(-1.0952605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2367661) q[0];
sx q[0];
rz(-2.5425826) q[0];
sx q[0];
rz(-2.616171) q[0];
x q[1];
rz(-0.11028408) q[2];
sx q[2];
rz(-0.37026893) q[2];
sx q[2];
rz(-0.25843378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.238823) q[1];
sx q[1];
rz(-2.2939957) q[1];
sx q[1];
rz(-0.1478424) q[1];
x q[2];
rz(0.59703897) q[3];
sx q[3];
rz(-1.9970915) q[3];
sx q[3];
rz(-2.7753914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1209391) q[2];
sx q[2];
rz(-1.8755269) q[2];
sx q[2];
rz(0.90547639) q[2];
rz(0.71546537) q[3];
sx q[3];
rz(-2.4162636) q[3];
sx q[3];
rz(-0.65413094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7249971) q[0];
sx q[0];
rz(-1.7323957) q[0];
sx q[0];
rz(0.61687627) q[0];
rz(1.1299409) q[1];
sx q[1];
rz(-1.2810974) q[1];
sx q[1];
rz(1.7249736) q[1];
rz(2.2466725) q[2];
sx q[2];
rz(-0.85325675) q[2];
sx q[2];
rz(0.85430145) q[2];
rz(2.1917584) q[3];
sx q[3];
rz(-2.8845308) q[3];
sx q[3];
rz(-1.844248) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
