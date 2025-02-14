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
rz(2.8994695) q[0];
sx q[0];
rz(-2.4162633) q[0];
sx q[0];
rz(-0.90670937) q[0];
rz(2.6746305) q[1];
sx q[1];
rz(-0.90277201) q[1];
sx q[1];
rz(-2.8977107) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4451673) q[0];
sx q[0];
rz(-1.1761242) q[0];
sx q[0];
rz(0.68035521) q[0];
rz(-pi) q[1];
rz(0.68274926) q[2];
sx q[2];
rz(-1.054316) q[2];
sx q[2];
rz(0.7841332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0816312) q[1];
sx q[1];
rz(-1.75358) q[1];
sx q[1];
rz(1.6441397) q[1];
rz(-pi) q[2];
rz(0.21190819) q[3];
sx q[3];
rz(-0.82900199) q[3];
sx q[3];
rz(-0.96860368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4563518) q[2];
sx q[2];
rz(-2.6754003) q[2];
sx q[2];
rz(2.1166128) q[2];
rz(3.0023365) q[3];
sx q[3];
rz(-1.1956297) q[3];
sx q[3];
rz(-0.16417424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3716632) q[0];
sx q[0];
rz(-2.0491845) q[0];
sx q[0];
rz(-1.9364233) q[0];
rz(1.6734164) q[1];
sx q[1];
rz(-1.877715) q[1];
sx q[1];
rz(-1.42043) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2989216) q[0];
sx q[0];
rz(-2.3793829) q[0];
sx q[0];
rz(1.803714) q[0];
x q[1];
rz(1.7718138) q[2];
sx q[2];
rz(-0.2699142) q[2];
sx q[2];
rz(2.6156099) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.08733315) q[1];
sx q[1];
rz(-1.9078322) q[1];
sx q[1];
rz(-2.9059306) q[1];
rz(2.0492433) q[3];
sx q[3];
rz(-0.88913267) q[3];
sx q[3];
rz(-2.9605856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18985192) q[2];
sx q[2];
rz(-2.9517089) q[2];
sx q[2];
rz(0.80113775) q[2];
rz(2.7028911) q[3];
sx q[3];
rz(-1.8935685) q[3];
sx q[3];
rz(-2.0285105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7640215) q[0];
sx q[0];
rz(-0.29733297) q[0];
sx q[0];
rz(-2.0024894) q[0];
rz(0.26690075) q[1];
sx q[1];
rz(-1.2511988) q[1];
sx q[1];
rz(-2.7762754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4802088) q[0];
sx q[0];
rz(-1.4719506) q[0];
sx q[0];
rz(2.6264725) q[0];
rz(2.2659698) q[2];
sx q[2];
rz(-0.77689028) q[2];
sx q[2];
rz(-0.50621835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2899785) q[1];
sx q[1];
rz(-1.9383926) q[1];
sx q[1];
rz(-0.062003597) q[1];
rz(0.86287873) q[3];
sx q[3];
rz(-2.2886057) q[3];
sx q[3];
rz(0.77117111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.9109362) q[2];
sx q[2];
rz(-2.4012884) q[2];
sx q[2];
rz(0.12348565) q[2];
rz(-0.89928818) q[3];
sx q[3];
rz(-1.6920009) q[3];
sx q[3];
rz(-1.9062769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3673636) q[0];
sx q[0];
rz(-1.6772567) q[0];
sx q[0];
rz(-0.88957077) q[0];
rz(-2.3956237) q[1];
sx q[1];
rz(-1.4624701) q[1];
sx q[1];
rz(1.3847345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9002514) q[0];
sx q[0];
rz(-1.5129733) q[0];
sx q[0];
rz(-3.1067294) q[0];
rz(-0.32082243) q[2];
sx q[2];
rz(-0.66705441) q[2];
sx q[2];
rz(-0.49122444) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3900323) q[1];
sx q[1];
rz(-2.1728325) q[1];
sx q[1];
rz(1.7979509) q[1];
x q[2];
rz(3.0880726) q[3];
sx q[3];
rz(-0.53456261) q[3];
sx q[3];
rz(-3.0381615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7437462) q[2];
sx q[2];
rz(-2.3036239) q[2];
sx q[2];
rz(-2.1261334) q[2];
rz(0.022653496) q[3];
sx q[3];
rz(-1.7706784) q[3];
sx q[3];
rz(-0.13815752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37462336) q[0];
sx q[0];
rz(-1.8648819) q[0];
sx q[0];
rz(-2.936777) q[0];
rz(0.21518937) q[1];
sx q[1];
rz(-0.22061017) q[1];
sx q[1];
rz(-0.87517103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3615205) q[0];
sx q[0];
rz(-0.339739) q[0];
sx q[0];
rz(1.1113412) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.740848) q[2];
sx q[2];
rz(-2.399246) q[2];
sx q[2];
rz(2.5352853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2339911) q[1];
sx q[1];
rz(-2.7669766) q[1];
sx q[1];
rz(-2.5020155) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93980333) q[3];
sx q[3];
rz(-1.4550337) q[3];
sx q[3];
rz(-0.94652688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.63336664) q[2];
sx q[2];
rz(-2.2534011) q[2];
sx q[2];
rz(1.4166895) q[2];
rz(-0.9238227) q[3];
sx q[3];
rz(-0.97771907) q[3];
sx q[3];
rz(-2.0098604) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1205587) q[0];
sx q[0];
rz(-0.73796213) q[0];
sx q[0];
rz(0.84838867) q[0];
rz(-1.5658762) q[1];
sx q[1];
rz(-1.8137685) q[1];
sx q[1];
rz(2.1956992) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28117546) q[0];
sx q[0];
rz(-2.3957154) q[0];
sx q[0];
rz(-0.4263878) q[0];
rz(1.7849493) q[2];
sx q[2];
rz(-2.3698167) q[2];
sx q[2];
rz(-2.0146148) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1887363) q[1];
sx q[1];
rz(-0.99299255) q[1];
sx q[1];
rz(-2.0835447) q[1];
rz(-pi) q[2];
rz(1.3729783) q[3];
sx q[3];
rz(-2.0536011) q[3];
sx q[3];
rz(-1.609926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9653603) q[2];
sx q[2];
rz(-1.9051899) q[2];
sx q[2];
rz(-1.7410834) q[2];
rz(1.4817574) q[3];
sx q[3];
rz(-1.7912495) q[3];
sx q[3];
rz(-2.9422133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4014715) q[0];
sx q[0];
rz(-0.49982247) q[0];
sx q[0];
rz(-1.9299782) q[0];
rz(0.47677332) q[1];
sx q[1];
rz(-1.8649273) q[1];
sx q[1];
rz(1.5062987) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9325628) q[0];
sx q[0];
rz(-1.8537921) q[0];
sx q[0];
rz(-2.1542633) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8938619) q[2];
sx q[2];
rz(-1.3328716) q[2];
sx q[2];
rz(2.0428773) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.41667485) q[1];
sx q[1];
rz(-0.70958269) q[1];
sx q[1];
rz(2.9516944) q[1];
rz(2.6166418) q[3];
sx q[3];
rz(-0.54407507) q[3];
sx q[3];
rz(1.8067738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8022884) q[2];
sx q[2];
rz(-1.6243287) q[2];
sx q[2];
rz(0.18292546) q[2];
rz(-2.4798992) q[3];
sx q[3];
rz(-1.9293834) q[3];
sx q[3];
rz(-2.4367512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6905007) q[0];
sx q[0];
rz(-1.4649614) q[0];
sx q[0];
rz(0.14725421) q[0];
rz(-2.9946949) q[1];
sx q[1];
rz(-2.2203827) q[1];
sx q[1];
rz(-1.9619092) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058181949) q[0];
sx q[0];
rz(-1.1549875) q[0];
sx q[0];
rz(1.2021628) q[0];
rz(2.0496164) q[2];
sx q[2];
rz(-1.1767595) q[2];
sx q[2];
rz(-1.0543329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50936401) q[1];
sx q[1];
rz(-1.6809789) q[1];
sx q[1];
rz(-2.5446507) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3875119) q[3];
sx q[3];
rz(-1.4907661) q[3];
sx q[3];
rz(-0.66129337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94072056) q[2];
sx q[2];
rz(-2.2042037) q[2];
sx q[2];
rz(0.76784039) q[2];
rz(2.614894) q[3];
sx q[3];
rz(-1.0517164) q[3];
sx q[3];
rz(0.55829486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.679477) q[0];
sx q[0];
rz(-0.18590346) q[0];
sx q[0];
rz(1.0900849) q[0];
rz(1.7915626) q[1];
sx q[1];
rz(-1.7410024) q[1];
sx q[1];
rz(0.42492351) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37575484) q[0];
sx q[0];
rz(-2.0415487) q[0];
sx q[0];
rz(-3.077917) q[0];
rz(-pi) q[1];
rz(0.64984729) q[2];
sx q[2];
rz(-1.51553) q[2];
sx q[2];
rz(-1.856232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7577014) q[1];
sx q[1];
rz(-0.5846068) q[1];
sx q[1];
rz(-0.65442185) q[1];
x q[2];
rz(0.70885371) q[3];
sx q[3];
rz(-1.9836622) q[3];
sx q[3];
rz(-2.9065437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2316042) q[2];
sx q[2];
rz(-0.51077545) q[2];
sx q[2];
rz(-2.6166022) q[2];
rz(-1.7867583) q[3];
sx q[3];
rz(-2.2460008) q[3];
sx q[3];
rz(1.7790599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43596426) q[0];
sx q[0];
rz(-0.66052496) q[0];
sx q[0];
rz(3.0349162) q[0];
rz(1.0542997) q[1];
sx q[1];
rz(-1.7091457) q[1];
sx q[1];
rz(-1.4073102) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5085691) q[0];
sx q[0];
rz(-0.79344024) q[0];
sx q[0];
rz(0.88280789) q[0];
x q[1];
rz(-0.65347341) q[2];
sx q[2];
rz(-1.1879041) q[2];
sx q[2];
rz(-3.058397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.52508721) q[1];
sx q[1];
rz(-2.8697851) q[1];
sx q[1];
rz(-2.1735783) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.484177) q[3];
sx q[3];
rz(-0.86674009) q[3];
sx q[3];
rz(0.92368607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4878896) q[2];
sx q[2];
rz(-2.6370878) q[2];
sx q[2];
rz(2.8686236) q[2];
rz(2.2702787) q[3];
sx q[3];
rz(-1.4600236) q[3];
sx q[3];
rz(-3.0070846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3684261) q[0];
sx q[0];
rz(-1.9482524) q[0];
sx q[0];
rz(0.85309749) q[0];
rz(0.38289616) q[1];
sx q[1];
rz(-1.2877512) q[1];
sx q[1];
rz(3.126694) q[1];
rz(0.56131526) q[2];
sx q[2];
rz(-1.9519162) q[2];
sx q[2];
rz(-1.219195) q[2];
rz(3.0440108) q[3];
sx q[3];
rz(-2.3392936) q[3];
sx q[3];
rz(0.38264075) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
