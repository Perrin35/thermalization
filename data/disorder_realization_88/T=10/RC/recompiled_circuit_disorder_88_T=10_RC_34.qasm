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
rz(0.091436401) q[1];
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
rz(0.63011516) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.66274553) q[1];
sx q[1];
rz(-1.6531684) q[1];
sx q[1];
rz(-1.3326416) q[1];
x q[2];
rz(-1.2535291) q[3];
sx q[3];
rz(-0.61875611) q[3];
sx q[3];
rz(1.3878824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4771007) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0098410957) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(-0.69357187) q[0];
rz(-2.0454848) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(0.19031659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41102558) q[0];
sx q[0];
rz(-1.1682296) q[0];
sx q[0];
rz(-2.9257724) q[0];
x q[1];
rz(1.9637738) q[2];
sx q[2];
rz(-1.0699546) q[2];
sx q[2];
rz(0.21265342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39365921) q[1];
sx q[1];
rz(-2.2927082) q[1];
sx q[1];
rz(0.02971239) q[1];
rz(0.5204366) q[3];
sx q[3];
rz(-0.79870975) q[3];
sx q[3];
rz(0.48898104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(1.2197536) q[2];
rz(-0.1427342) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74293566) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(0.29207686) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0322198) q[0];
sx q[0];
rz(-1.4782227) q[0];
sx q[0];
rz(-2.5412482) q[0];
x q[1];
rz(1.5065932) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(-1.1484255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51860147) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(1.4856505) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1850584) q[3];
sx q[3];
rz(-2.4534561) q[3];
sx q[3];
rz(-3.0582173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-2.0911066) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7320025) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(-0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.2971372) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521833) q[0];
sx q[0];
rz(-1.4404738) q[0];
sx q[0];
rz(-0.043686314) q[0];
rz(-pi) q[1];
rz(-1.2816216) q[2];
sx q[2];
rz(-1.4099858) q[2];
sx q[2];
rz(-0.5493872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2343826) q[1];
sx q[1];
rz(-2.130059) q[1];
sx q[1];
rz(0.84749605) q[1];
rz(1.6468871) q[3];
sx q[3];
rz(-0.18008672) q[3];
sx q[3];
rz(2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90594784) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(2.2223991) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.250524) q[0];
sx q[0];
rz(-0.91618012) q[0];
sx q[0];
rz(-3.0380681) q[0];
rz(-pi) q[1];
rz(1.2582785) q[2];
sx q[2];
rz(-1.7973571) q[2];
sx q[2];
rz(0.3414008) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1713472) q[1];
sx q[1];
rz(-0.73927021) q[1];
sx q[1];
rz(1.7319748) q[1];
rz(-pi) q[2];
rz(-0.77107314) q[3];
sx q[3];
rz(-1.8043451) q[3];
sx q[3];
rz(1.1119103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(0.041794725) q[2];
rz(-3.0801008) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(2.1110995) q[0];
rz(2.4018535) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(2.5700263) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58493462) q[0];
sx q[0];
rz(-0.99781636) q[0];
sx q[0];
rz(-2.1893188) q[0];
x q[1];
rz(0.12708725) q[2];
sx q[2];
rz(-0.90887827) q[2];
sx q[2];
rz(1.6866637) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73733444) q[1];
sx q[1];
rz(-1.8109651) q[1];
sx q[1];
rz(1.1024229) q[1];
rz(2.3011175) q[3];
sx q[3];
rz(-0.90208902) q[3];
sx q[3];
rz(0.83276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.968154) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(-2.5773933) q[2];
rz(-3.0155904) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(-2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(2.4760822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3367046) q[0];
sx q[0];
rz(-2.8993653) q[0];
sx q[0];
rz(-1.5664943) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.268157) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(2.9261677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1674541) q[1];
sx q[1];
rz(-2.0270837) q[1];
sx q[1];
rz(-2.8970701) q[1];
rz(1.8129559) q[3];
sx q[3];
rz(-1.7112964) q[3];
sx q[3];
rz(1.5618009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9843288) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(1.4228014) q[2];
rz(3.026399) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(-2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049906235) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(-2.514839) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-2.802882) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9324547) q[0];
sx q[0];
rz(-1.4203686) q[0];
sx q[0];
rz(-1.223279) q[0];
rz(1.8024826) q[2];
sx q[2];
rz(-2.0121687) q[2];
sx q[2];
rz(2.1053134) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.415886) q[1];
sx q[1];
rz(-2.6542414) q[1];
sx q[1];
rz(1.332167) q[1];
x q[2];
rz(-2.268928) q[3];
sx q[3];
rz(-2.7250395) q[3];
sx q[3];
rz(2.0779028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64615858) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(0.70043606) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(0.65834808) q[0];
rz(2.530653) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(0.13959612) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0567719) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(1.8875185) q[0];
rz(-pi) q[1];
rz(-2.326968) q[2];
sx q[2];
rz(-1.5280208) q[2];
sx q[2];
rz(0.2864366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4715251) q[1];
sx q[1];
rz(-1.0499665) q[1];
sx q[1];
rz(-1.6478369) q[1];
x q[2];
rz(-0.052501909) q[3];
sx q[3];
rz(-0.92696654) q[3];
sx q[3];
rz(0.98464636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6167986) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(-0.33111462) q[2];
rz(0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-3.0537135) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(-0.18558003) q[0];
rz(2.045385) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(-1.6569482) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164924) q[0];
sx q[0];
rz(-2.6161368) q[0];
sx q[0];
rz(1.0111965) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9160203) q[2];
sx q[2];
rz(-1.8403056) q[2];
sx q[2];
rz(-0.05664209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2828196) q[1];
sx q[1];
rz(-0.38837896) q[1];
sx q[1];
rz(-1.8358843) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5649892) q[3];
sx q[3];
rz(-0.70492893) q[3];
sx q[3];
rz(2.4095636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93402702) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(2.5893842) q[2];
rz(2.3637555) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375297) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(2.8616703) q[2];
sx q[2];
rz(-1.7114077) q[2];
sx q[2];
rz(-0.68243295) q[2];
rz(-0.68708146) q[3];
sx q[3];
rz(-1.0071181) q[3];
sx q[3];
rz(-1.3345171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];