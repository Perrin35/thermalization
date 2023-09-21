OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(-0.29247984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777893) q[0];
sx q[0];
rz(-1.6262494) q[0];
sx q[0];
rz(1.1763014) q[0];
rz(1.6220665) q[2];
sx q[2];
rz(-2.3495557) q[2];
sx q[2];
rz(-0.37026065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3136183) q[1];
sx q[1];
rz(-0.59578124) q[1];
sx q[1];
rz(-1.4890563) q[1];
rz(-2.0446159) q[3];
sx q[3];
rz(-1.5690656) q[3];
sx q[3];
rz(0.8918106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.73074377) q[2];
sx q[2];
rz(-1.0935254) q[2];
sx q[2];
rz(-0.4494108) q[2];
rz(-0.64569008) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(-1.2908363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9280424) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(0.74203062) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.7785633) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33107685) q[0];
sx q[0];
rz(-1.1805981) q[0];
sx q[0];
rz(-0.66731989) q[0];
x q[1];
rz(1.7861373) q[2];
sx q[2];
rz(-1.1110348) q[2];
sx q[2];
rz(-0.41972566) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.084761707) q[1];
sx q[1];
rz(-0.34013224) q[1];
sx q[1];
rz(-0.057573307) q[1];
rz(0.19248776) q[3];
sx q[3];
rz(-1.3093595) q[3];
sx q[3];
rz(-2.1646433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8403975) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.8691241) q[2];
rz(-0.33723351) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(-0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0740046) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(2.8787676) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(-1.9062818) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6567898) q[0];
sx q[0];
rz(-3.1233388) q[0];
sx q[0];
rz(-0.59469964) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.026651816) q[2];
sx q[2];
rz(-1.3304552) q[2];
sx q[2];
rz(1.9453366) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1097818) q[1];
sx q[1];
rz(-0.33966741) q[1];
sx q[1];
rz(-2.8207645) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2658306) q[3];
sx q[3];
rz(-1.8417629) q[3];
sx q[3];
rz(2.392644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-0.18033218) q[2];
rz(0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(-2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76698774) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(-0.71907991) q[0];
rz(-2.6285697) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(-1.0346574) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35095222) q[0];
sx q[0];
rz(-1.5119023) q[0];
sx q[0];
rz(1.590593) q[0];
x q[1];
rz(-1.9070508) q[2];
sx q[2];
rz(-2.9985235) q[2];
sx q[2];
rz(-2.2640995) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0532916) q[1];
sx q[1];
rz(-0.75690714) q[1];
sx q[1];
rz(0.41815586) q[1];
rz(-pi) q[2];
rz(-1.7454809) q[3];
sx q[3];
rz(-1.1038728) q[3];
sx q[3];
rz(-0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2085312) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(2.9053524) q[2];
rz(1.9832206) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(-0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(-2.2205655) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-2.5193118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.92684) q[0];
sx q[0];
rz(-1.8536957) q[0];
sx q[0];
rz(-0.19006417) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4786517) q[2];
sx q[2];
rz(-2.2170057) q[2];
sx q[2];
rz(0.19942936) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29187782) q[1];
sx q[1];
rz(-1.4018702) q[1];
sx q[1];
rz(0.090976322) q[1];
rz(-pi) q[2];
rz(-2.4683687) q[3];
sx q[3];
rz(-2.554318) q[3];
sx q[3];
rz(-0.67929635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2174012) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(-3.1203111) q[2];
rz(1.6993258) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845881) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(0.044629991) q[0];
rz(2.1394829) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1761985) q[0];
sx q[0];
rz(-1.5571496) q[0];
sx q[0];
rz(1.3854881) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7897286) q[2];
sx q[2];
rz(-0.5133252) q[2];
sx q[2];
rz(2.9190612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1019563) q[1];
sx q[1];
rz(-2.4974681) q[1];
sx q[1];
rz(-1.0108175) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4061635) q[3];
sx q[3];
rz(-2.2693686) q[3];
sx q[3];
rz(-0.64100953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78952152) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(0.96674031) q[2];
rz(3.1294075) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(-0.10908443) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352585) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(0.99408856) q[0];
rz(-3.0864691) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(0.73928839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3294322) q[0];
sx q[0];
rz(-1.433916) q[0];
sx q[0];
rz(-0.095892266) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57619608) q[2];
sx q[2];
rz(-1.155876) q[2];
sx q[2];
rz(0.31838271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97923219) q[1];
sx q[1];
rz(-1.7602966) q[1];
sx q[1];
rz(2.6081671) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39412734) q[3];
sx q[3];
rz(-2.4836802) q[3];
sx q[3];
rz(3.1042299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5560975) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(-2.5778256) q[2];
rz(0.24027763) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(-0.58404303) q[0];
rz(-0.42075992) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(-2.8631794) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95094341) q[0];
sx q[0];
rz(-1.4454495) q[0];
sx q[0];
rz(2.1423253) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39324795) q[2];
sx q[2];
rz(-1.3967782) q[2];
sx q[2];
rz(-0.72088036) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.091374) q[1];
sx q[1];
rz(-0.41793567) q[1];
sx q[1];
rz(2.1347079) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9592459) q[3];
sx q[3];
rz(-0.58245814) q[3];
sx q[3];
rz(0.1633446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86928308) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(-2.753624) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(-2.8665682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824654) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(-0.7318837) q[0];
rz(-0.16648079) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(-0.073721185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1151704) q[0];
sx q[0];
rz(-2.0048884) q[0];
sx q[0];
rz(0.11154453) q[0];
rz(-pi) q[1];
rz(-1.45103) q[2];
sx q[2];
rz(-1.6511859) q[2];
sx q[2];
rz(-2.5188841) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32578711) q[1];
sx q[1];
rz(-2.1539638) q[1];
sx q[1];
rz(-1.6156242) q[1];
rz(-pi) q[2];
rz(-0.57314408) q[3];
sx q[3];
rz(-1.4775652) q[3];
sx q[3];
rz(1.4179109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5471197) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(2.7129042) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(1.2111506) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73228943) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(-2.0987341) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(1.7932549) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74603426) q[0];
sx q[0];
rz(-1.4575882) q[0];
sx q[0];
rz(-0.025870196) q[0];
rz(-pi) q[1];
rz(-1.4820319) q[2];
sx q[2];
rz(-2.4030493) q[2];
sx q[2];
rz(0.86631394) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18689449) q[1];
sx q[1];
rz(-2.8098626) q[1];
sx q[1];
rz(-1.3518672) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24095778) q[3];
sx q[3];
rz(-1.7562508) q[3];
sx q[3];
rz(1.0518187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(-0.094816118) q[2];
rz(-1.9684277) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6288347) q[0];
sx q[0];
rz(-2.6678968) q[0];
sx q[0];
rz(-2.0976023) q[0];
rz(-1.5402773) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(3.0586254) q[2];
sx q[2];
rz(-1.3583747) q[2];
sx q[2];
rz(2.8521392) q[2];
rz(0.81090609) q[3];
sx q[3];
rz(-1.1699642) q[3];
sx q[3];
rz(-1.5311833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];