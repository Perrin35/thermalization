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
rz(1.1458122) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(2.8491128) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4300633) q[0];
sx q[0];
rz(-1.9646514) q[0];
sx q[0];
rz(-3.0815365) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3621759) q[2];
sx q[2];
rz(-1.6072818) q[2];
sx q[2];
rz(1.905029) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3136183) q[1];
sx q[1];
rz(-2.5458114) q[1];
sx q[1];
rz(-1.6525364) q[1];
x q[2];
rz(3.1396477) q[3];
sx q[3];
rz(-1.0969775) q[3];
sx q[3];
rz(0.67809826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4108489) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(-0.4494108) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(1.2908363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9280424) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(2.399562) q[0];
rz(-1.4713326) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.3630294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5308967) q[0];
sx q[0];
rz(-0.96141059) q[0];
sx q[0];
rz(-1.0884398) q[0];
x q[1];
rz(2.6724733) q[2];
sx q[2];
rz(-1.763478) q[2];
sx q[2];
rz(-1.0543146) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7098397) q[1];
sx q[1];
rz(-1.5515986) q[1];
sx q[1];
rz(-0.33961105) q[1];
rz(-pi) q[2];
rz(0.19248776) q[3];
sx q[3];
rz(-1.3093595) q[3];
sx q[3];
rz(-2.1646433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8403975) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(-1.2724686) q[2];
rz(-0.33723351) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(0.2628251) q[0];
rz(0.35573959) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(-1.9062818) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6806157) q[0];
sx q[0];
rz(-1.5605698) q[0];
sx q[0];
rz(3.1264722) q[0];
rz(-pi) q[1];
rz(-3.1149408) q[2];
sx q[2];
rz(-1.8111374) q[2];
sx q[2];
rz(1.9453366) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.77093) q[1];
sx q[1];
rz(-1.2491033) q[1];
sx q[1];
rz(1.6817723) q[1];
rz(0.87576207) q[3];
sx q[3];
rz(-1.2998298) q[3];
sx q[3];
rz(-2.392644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(2.9612605) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.9790117) q[3];
sx q[3];
rz(0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3746049) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(-0.71907991) q[0];
rz(2.6285697) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(-1.0346574) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7906404) q[0];
sx q[0];
rz(-1.6296903) q[0];
sx q[0];
rz(-1.590593) q[0];
rz(1.7059533) q[2];
sx q[2];
rz(-1.617859) q[2];
sx q[2];
rz(1.0263718) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82995854) q[1];
sx q[1];
rz(-1.8533851) q[1];
sx q[1];
rz(-0.71210536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8095874) q[3];
sx q[3];
rz(-0.49626353) q[3];
sx q[3];
rz(0.23976025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93306142) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(2.9053524) q[2];
rz(1.9832206) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(2.2896144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(-2.239256) q[0];
rz(-0.92102712) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-0.62228084) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21475269) q[0];
sx q[0];
rz(-1.2878969) q[0];
sx q[0];
rz(-0.19006417) q[0];
rz(-pi) q[1];
rz(1.4786517) q[2];
sx q[2];
rz(-0.92458692) q[2];
sx q[2];
rz(-2.9421633) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8497148) q[1];
sx q[1];
rz(-1.4018702) q[1];
sx q[1];
rz(0.090976322) q[1];
x q[2];
rz(2.6617674) q[3];
sx q[3];
rz(-1.92355) q[3];
sx q[3];
rz(1.6640116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2174012) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(-3.1203111) q[2];
rz(-1.6993258) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(2.2309979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845881) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(-3.0969627) q[0];
rz(-1.0021098) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(3.1071641) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653942) q[0];
sx q[0];
rz(-1.5571496) q[0];
sx q[0];
rz(-1.7561046) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0197633) q[2];
sx q[2];
rz(-1.0708772) q[2];
sx q[2];
rz(-3.1140529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1019563) q[1];
sx q[1];
rz(-2.4974681) q[1];
sx q[1];
rz(2.1307751) q[1];
rz(-0.19272007) q[3];
sx q[3];
rz(-0.7145213) q[3];
sx q[3];
rz(-0.89380985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3520711) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(-0.96674031) q[2];
rz(-0.012185193) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(-3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.1352585) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(-0.99408856) q[0];
rz(-3.0864691) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(-2.4023043) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2544884) q[0];
sx q[0];
rz(-1.4758037) q[0];
sx q[0];
rz(1.4332921) q[0];
x q[1];
rz(2.5653966) q[2];
sx q[2];
rz(-1.9857166) q[2];
sx q[2];
rz(-0.31838271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48078295) q[1];
sx q[1];
rz(-1.0479095) q[1];
sx q[1];
rz(-1.3516264) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8592632) q[3];
sx q[3];
rz(-0.97086421) q[3];
sx q[3];
rz(-0.44655061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5560975) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(-0.24027763) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(-2.5575496) q[0];
rz(2.7208327) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(-0.27841321) q[1];
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
x q[1];
rz(2.7114696) q[2];
sx q[2];
rz(-2.7134036) q[2];
sx q[2];
rz(2.6870514) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.091374) q[1];
sx q[1];
rz(-0.41793567) q[1];
sx q[1];
rz(-1.0068847) q[1];
rz(-pi) q[2];
rz(1.0233381) q[3];
sx q[3];
rz(-1.7806782) q[3];
sx q[3];
rz(-1.7367712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2723096) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(-0.38796866) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(2.8665682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85912722) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(-2.4097089) q[0];
rz(0.16648079) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(3.0678715) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28669824) q[0];
sx q[0];
rz(-2.694283) q[0];
sx q[0];
rz(1.8064503) q[0];
x q[1];
rz(-3.0606255) q[2];
sx q[2];
rz(-1.4514187) q[2];
sx q[2];
rz(-0.93842426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2697061) q[1];
sx q[1];
rz(-1.5333813) q[1];
sx q[1];
rz(0.58362959) q[1];
rz(1.681626) q[3];
sx q[3];
rz(-1.0004527) q[3];
sx q[3];
rz(-0.092872083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5471197) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(2.7129042) q[2];
rz(1.3171014) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4093032) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(-2.0987341) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(1.3483378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52085984) q[0];
sx q[0];
rz(-0.11611406) q[0];
sx q[0];
rz(-1.7945047) q[0];
rz(-1.6595608) q[2];
sx q[2];
rz(-0.73854337) q[2];
sx q[2];
rz(-2.2752787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1765602) q[1];
sx q[1];
rz(-1.6415879) q[1];
sx q[1];
rz(-1.895158) q[1];
x q[2];
rz(2.9006349) q[3];
sx q[3];
rz(-1.7562508) q[3];
sx q[3];
rz(-2.0897739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80031359) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(3.0467765) q[2];
rz(-1.9684277) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(-0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6288347) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(1.6013153) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(1.2039456) q[2];
sx q[2];
rz(-2.9137712) q[2];
sx q[2];
rz(-3.055345) q[2];
rz(-0.81090609) q[3];
sx q[3];
rz(-1.9716284) q[3];
sx q[3];
rz(1.6104094) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];