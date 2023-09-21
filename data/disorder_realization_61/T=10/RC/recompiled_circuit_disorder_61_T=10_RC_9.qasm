OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(1.1336741) q[0];
sx q[0];
rz(7.8757474) q[0];
rz(1.6917317) q[1];
sx q[1];
rz(-0.65728846) q[1];
sx q[1];
rz(-2.5974098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.777594) q[0];
sx q[0];
rz(-1.5033804) q[0];
sx q[0];
rz(1.534034) q[0];
x q[1];
rz(1.5266225) q[2];
sx q[2];
rz(-1.4695115) q[2];
sx q[2];
rz(-0.12461187) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4197598) q[1];
sx q[1];
rz(-2.3082323) q[1];
sx q[1];
rz(1.014773) q[1];
x q[2];
rz(1.5531494) q[3];
sx q[3];
rz(-1.5421765) q[3];
sx q[3];
rz(2.0538405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.071775285) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(-1.7791746) q[2];
rz(-0.028256265) q[3];
sx q[3];
rz(-1.7621721) q[3];
sx q[3];
rz(2.3513667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966999) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(1.171296) q[0];
rz(-2.9303739) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(-1.7374977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1275741) q[0];
sx q[0];
rz(-1.1496135) q[0];
sx q[0];
rz(-0.76571) q[0];
rz(-pi) q[1];
rz(2.6589083) q[2];
sx q[2];
rz(-1.7182554) q[2];
sx q[2];
rz(0.54373103) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1120226) q[1];
sx q[1];
rz(-2.8575172) q[1];
sx q[1];
rz(0.082138852) q[1];
x q[2];
rz(-2.6614463) q[3];
sx q[3];
rz(-1.4501791) q[3];
sx q[3];
rz(3.0062825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8319548) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(0.94397604) q[2];
rz(-0.55654636) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(-1.8054307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495162) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(-1.4235494) q[0];
rz(0.81958333) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(-2.5779285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(1.4865727) q[0];
rz(-pi) q[1];
rz(-2.5416449) q[2];
sx q[2];
rz(-1.7581345) q[2];
sx q[2];
rz(-1.7768605) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3445661) q[1];
sx q[1];
rz(-2.1246506) q[1];
sx q[1];
rz(-0.38573854) q[1];
x q[2];
rz(2.5936801) q[3];
sx q[3];
rz(-0.69763819) q[3];
sx q[3];
rz(-1.9426949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6391969) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(1.0602661) q[2];
rz(3.0660196) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(-0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183427) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(-0.19317214) q[0];
rz(-1.5441783) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(-2.4838122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7155899) q[0];
sx q[0];
rz(-1.3960739) q[0];
sx q[0];
rz(1.6013553) q[0];
rz(-pi) q[1];
rz(-2.0006318) q[2];
sx q[2];
rz(-1.2865215) q[2];
sx q[2];
rz(-2.5179799) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2967589) q[1];
sx q[1];
rz(-1.3169603) q[1];
sx q[1];
rz(-1.0253419) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2781906) q[3];
sx q[3];
rz(-1.3232104) q[3];
sx q[3];
rz(2.2729371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0956991) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(1.0470541) q[2];
rz(-1.0632769) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(-1.261238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.7970153) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(2.1176594) q[0];
rz(2.3539885) q[1];
sx q[1];
rz(-1.5757631) q[1];
sx q[1];
rz(-2.5147298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9526564) q[0];
sx q[0];
rz(-3.1075826) q[0];
sx q[0];
rz(2.0859499) q[0];
x q[1];
rz(0.31693964) q[2];
sx q[2];
rz(-1.8570515) q[2];
sx q[2];
rz(2.0010819) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31506854) q[1];
sx q[1];
rz(-1.8693923) q[1];
sx q[1];
rz(1.8930757) q[1];
rz(-pi) q[2];
rz(-1.9634788) q[3];
sx q[3];
rz(-1.6851478) q[3];
sx q[3];
rz(-1.3254195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91288599) q[2];
sx q[2];
rz(-1.2330981) q[2];
sx q[2];
rz(1.0181001) q[2];
rz(-0.61156887) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(-2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927032) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(1.1992136) q[0];
rz(-1.9723643) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(2.8170524) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51844653) q[0];
sx q[0];
rz(-1.070797) q[0];
sx q[0];
rz(-0.48796939) q[0];
rz(-pi) q[1];
rz(2.4915699) q[2];
sx q[2];
rz(-0.92009244) q[2];
sx q[2];
rz(-2.0536154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2734087) q[1];
sx q[1];
rz(-1.5460099) q[1];
sx q[1];
rz(1.0367869) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8683897) q[3];
sx q[3];
rz(-1.9801567) q[3];
sx q[3];
rz(0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.2754053) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-2.349699) q[3];
sx q[3];
rz(1.5462497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642898) q[0];
sx q[0];
rz(-1.807656) q[0];
sx q[0];
rz(1.7328847) q[0];
rz(2.6858221) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(1.9394402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31583187) q[0];
sx q[0];
rz(-1.1734661) q[0];
sx q[0];
rz(-0.8855008) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7289594) q[2];
sx q[2];
rz(-0.79421439) q[2];
sx q[2];
rz(-2.7318294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68820885) q[1];
sx q[1];
rz(-1.6628478) q[1];
sx q[1];
rz(1.6545014) q[1];
x q[2];
rz(-2.8771411) q[3];
sx q[3];
rz(-1.7489986) q[3];
sx q[3];
rz(1.7314535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(0.84623519) q[2];
rz(2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071035944) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(-0.69560266) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(1.5323458) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42513645) q[0];
sx q[0];
rz(-1.813659) q[0];
sx q[0];
rz(-0.48432414) q[0];
rz(2.0267018) q[2];
sx q[2];
rz(-1.1258944) q[2];
sx q[2];
rz(-1.0783431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6032226) q[1];
sx q[1];
rz(-0.45062989) q[1];
sx q[1];
rz(1.5249114) q[1];
rz(-pi) q[2];
rz(-1.6892151) q[3];
sx q[3];
rz(-0.97202557) q[3];
sx q[3];
rz(2.1042202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(-2.2873986) q[2];
rz(1.2285852) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5478741) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(-2.4771931) q[0];
rz(1.9539072) q[1];
sx q[1];
rz(-2.4408051) q[1];
sx q[1];
rz(1.2845576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2537456) q[0];
sx q[0];
rz(-2.3652024) q[0];
sx q[0];
rz(1.1573769) q[0];
rz(-pi) q[1];
rz(1.9509435) q[2];
sx q[2];
rz(-1.322762) q[2];
sx q[2];
rz(1.5284577) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38332332) q[1];
sx q[1];
rz(-2.3136534) q[1];
sx q[1];
rz(0.62288706) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5434389) q[3];
sx q[3];
rz(-1.9606855) q[3];
sx q[3];
rz(-1.1297806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61218843) q[2];
sx q[2];
rz(-0.89451423) q[2];
sx q[2];
rz(-2.5908296) q[2];
rz(2.4380056) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4187014) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(0.044145949) q[0];
rz(-1.6126397) q[1];
sx q[1];
rz(-1.9020558) q[1];
sx q[1];
rz(-0.77967656) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75694376) q[0];
sx q[0];
rz(-2.1063519) q[0];
sx q[0];
rz(1.4063565) q[0];
x q[1];
rz(-0.81515628) q[2];
sx q[2];
rz(-1.2218468) q[2];
sx q[2];
rz(2.7100035) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70982198) q[1];
sx q[1];
rz(-0.55574544) q[1];
sx q[1];
rz(2.7906448) q[1];
rz(-2.6727647) q[3];
sx q[3];
rz(-2.3033597) q[3];
sx q[3];
rz(1.9459141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6471275) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(1.5987827) q[2];
rz(-0.70458448) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.71173944) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(-1.343887) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(1.1336397) q[2];
sx q[2];
rz(-1.2883678) q[2];
sx q[2];
rz(-1.012445) q[2];
rz(-1.703249) q[3];
sx q[3];
rz(-0.9369029) q[3];
sx q[3];
rz(2.423219) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
