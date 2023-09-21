OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(4.3972754) q[0];
sx q[0];
rz(9.7527405) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(-3.0501563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56337315) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(-2.5114775) q[0];
rz(-pi) q[1];
rz(0.25912164) q[2];
sx q[2];
rz(-1.4403575) q[2];
sx q[2];
rz(0.63333095) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5603197) q[1];
sx q[1];
rz(-2.8898507) q[1];
sx q[1];
rz(1.2341577) q[1];
rz(-pi) q[2];
rz(-1.2535291) q[3];
sx q[3];
rz(-0.61875611) q[3];
sx q[3];
rz(-1.7537102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(-1.1260024) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(-0.69357187) q[0];
rz(-1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(0.19031659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099146518) q[0];
sx q[0];
rz(-2.6876246) q[0];
sx q[0];
rz(-1.1048261) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53493494) q[2];
sx q[2];
rz(-1.2282279) q[2];
sx q[2];
rz(1.5869706) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.39365921) q[1];
sx q[1];
rz(-2.2927082) q[1];
sx q[1];
rz(0.02971239) q[1];
x q[2];
rz(-0.5204366) q[3];
sx q[3];
rz(-2.3428829) q[3];
sx q[3];
rz(-2.6526116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50743121) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(-0.1427342) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(-2.232961) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74293566) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(2.8495158) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(1.8935727) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0322198) q[0];
sx q[0];
rz(-1.4782227) q[0];
sx q[0];
rz(-0.60034445) q[0];
x q[1];
rz(-2.9343611) q[2];
sx q[2];
rz(-1.5079632) q[2];
sx q[2];
rz(-0.40916967) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51860147) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(1.4856505) q[1];
x q[2];
rz(-2.2218024) q[3];
sx q[3];
rz(-1.3295104) q[3];
sx q[3];
rz(1.9581219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(1.9281663) q[2];
rz(-2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40959013) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(1.2971372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064917795) q[0];
sx q[0];
rz(-0.13741048) q[0];
sx q[0];
rz(1.2491559) q[0];
rz(0.16764955) q[2];
sx q[2];
rz(-1.285458) q[2];
sx q[2];
rz(-0.97380762) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2343826) q[1];
sx q[1];
rz(-1.0115336) q[1];
sx q[1];
rz(2.2940966) q[1];
x q[2];
rz(1.6468871) q[3];
sx q[3];
rz(-0.18008672) q[3];
sx q[3];
rz(-1.0877346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2356448) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(2.2223991) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(-1.2478158) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(1.8163266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(-1.7153046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.250524) q[0];
sx q[0];
rz(-0.91618012) q[0];
sx q[0];
rz(-3.0380681) q[0];
x q[1];
rz(0.9274474) q[2];
sx q[2];
rz(-2.7577835) q[2];
sx q[2];
rz(2.5197033) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95477415) q[1];
sx q[1];
rz(-2.2983119) q[1];
sx q[1];
rz(2.9963042) q[1];
rz(-1.8910847) q[3];
sx q[3];
rz(-0.8257782) q[3];
sx q[3];
rz(2.461493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(-3.0997979) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(-2.7048236) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3549266) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(-1.0304931) q[0];
rz(2.4018535) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(-0.57156634) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067236) q[0];
sx q[0];
rz(-2.3248701) q[0];
sx q[0];
rz(2.4094894) q[0];
rz(-1.4095441) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(1.6598998) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4042582) q[1];
sx q[1];
rz(-1.3306276) q[1];
sx q[1];
rz(2.0391697) q[1];
rz(0.84047517) q[3];
sx q[3];
rz(-0.90208902) q[3];
sx q[3];
rz(2.3088282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.968154) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(-0.56419939) q[2];
rz(3.0155904) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(2.8469767) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903704) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(0.5258711) q[1];
sx q[1];
rz(-0.41627517) q[1];
sx q[1];
rz(0.66551048) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80488801) q[0];
sx q[0];
rz(-0.24222736) q[0];
sx q[0];
rz(-1.5664943) q[0];
x q[1];
rz(2.268157) q[2];
sx q[2];
rz(-1.6849815) q[2];
sx q[2];
rz(2.9261677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9741386) q[1];
sx q[1];
rz(-1.1145089) q[1];
sx q[1];
rz(0.2445226) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0379167) q[3];
sx q[3];
rz(-2.8623192) q[3];
sx q[3];
rz(-2.6168407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9843288) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(-3.026399) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049906235) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(-2.514839) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(-0.33871067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3876462) q[0];
sx q[0];
rz(-0.37746143) q[0];
sx q[0];
rz(-1.1520555) q[0];
x q[1];
rz(-0.4523925) q[2];
sx q[2];
rz(-2.6466742) q[2];
sx q[2];
rz(2.6099043) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.056811) q[1];
sx q[1];
rz(-1.6817131) q[1];
sx q[1];
rz(-2.0463498) q[1];
rz(-1.2440153) q[3];
sx q[3];
rz(-1.3076926) q[3];
sx q[3];
rz(2.9941878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4954341) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(-0.70043606) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(-0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-2.4832446) q[0];
rz(2.530653) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(-3.0019965) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7394373) q[0];
sx q[0];
rz(-1.6318775) q[0];
sx q[0];
rz(-0.020045965) q[0];
rz(3.082824) q[2];
sx q[2];
rz(-0.81548703) q[2];
sx q[2];
rz(-1.8975443) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86233625) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(-0.52211296) q[1];
x q[2];
rz(-3.0890907) q[3];
sx q[3];
rz(-2.2146261) q[3];
sx q[3];
rz(0.98464636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6167986) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(2.810478) q[2];
rz(-2.3838499) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-3.0537135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963294) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(-2.9560126) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.6569482) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0508704) q[0];
sx q[0];
rz(-1.301268) q[0];
sx q[0];
rz(-1.1140633) q[0];
rz(0.22557232) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(-3.0849506) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46586043) q[1];
sx q[1];
rz(-1.4714186) q[1];
sx q[1];
rz(-1.194721) q[1];
rz(-2.2757169) q[3];
sx q[3];
rz(-1.5670334) q[3];
sx q[3];
rz(-0.83434425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2075656) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(2.5893842) q[2];
rz(0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8778397) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(0.18763018) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(1.7170231) q[2];
sx q[2];
rz(-1.8478827) q[2];
sx q[2];
rz(-2.293496) q[2];
rz(-2.256176) q[3];
sx q[3];
rz(-2.1366742) q[3];
sx q[3];
rz(-2.4921806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];