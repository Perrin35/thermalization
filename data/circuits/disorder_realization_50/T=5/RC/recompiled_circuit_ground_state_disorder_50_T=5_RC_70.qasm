OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.1640373) q[0];
sx q[0];
rz(-1.9174175) q[0];
sx q[0];
rz(-1.2807711) q[0];
rz(-0.7723074) q[1];
sx q[1];
rz(-0.63900715) q[1];
sx q[1];
rz(0.26689902) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3506265) q[0];
sx q[0];
rz(-1.6668975) q[0];
sx q[0];
rz(3.0477357) q[0];
x q[1];
rz(0.25702567) q[2];
sx q[2];
rz(-0.89867175) q[2];
sx q[2];
rz(-0.66115236) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81921009) q[1];
sx q[1];
rz(-1.9967845) q[1];
sx q[1];
rz(0.8041348) q[1];
rz(-pi) q[2];
rz(2.9270615) q[3];
sx q[3];
rz(-0.45154587) q[3];
sx q[3];
rz(-2.8928234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.018517) q[2];
sx q[2];
rz(-2.1776431) q[2];
sx q[2];
rz(-0.71259552) q[2];
rz(0.16768843) q[3];
sx q[3];
rz(-0.84961397) q[3];
sx q[3];
rz(2.6970421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0406822) q[0];
sx q[0];
rz(-1.3586783) q[0];
sx q[0];
rz(-1.6810625) q[0];
rz(-1.679861) q[1];
sx q[1];
rz(-2.0748383) q[1];
sx q[1];
rz(1.1163968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8212962) q[0];
sx q[0];
rz(-1.5954994) q[0];
sx q[0];
rz(1.5276093) q[0];
x q[1];
rz(2.9903611) q[2];
sx q[2];
rz(-1.0200715) q[2];
sx q[2];
rz(3.0557003) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9560686) q[1];
sx q[1];
rz(-1.9466234) q[1];
sx q[1];
rz(1.8517428) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7313569) q[3];
sx q[3];
rz(-1.6950775) q[3];
sx q[3];
rz(1.4924442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1560893) q[2];
sx q[2];
rz(-2.6458793) q[2];
sx q[2];
rz(2.8893341) q[2];
rz(-1.0559399) q[3];
sx q[3];
rz(-2.2838433) q[3];
sx q[3];
rz(1.0004388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5591739) q[0];
sx q[0];
rz(-0.90455872) q[0];
sx q[0];
rz(-0.82758033) q[0];
rz(0.52455348) q[1];
sx q[1];
rz(-1.8962212) q[1];
sx q[1];
rz(0.96447271) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9081952) q[0];
sx q[0];
rz(-0.99401532) q[0];
sx q[0];
rz(-0.080349313) q[0];
x q[1];
rz(1.1954514) q[2];
sx q[2];
rz(-1.6520733) q[2];
sx q[2];
rz(-0.92513212) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2169184) q[1];
sx q[1];
rz(-2.5092832) q[1];
sx q[1];
rz(3.0576474) q[1];
rz(-pi) q[2];
rz(0.44154756) q[3];
sx q[3];
rz(-1.5184622) q[3];
sx q[3];
rz(-0.56091058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39782897) q[2];
sx q[2];
rz(-0.22481329) q[2];
sx q[2];
rz(1.1986097) q[2];
rz(-1.4926636) q[3];
sx q[3];
rz(-2.1677833) q[3];
sx q[3];
rz(-0.57884136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8156133) q[0];
sx q[0];
rz(-2.9878243) q[0];
sx q[0];
rz(1.9258668) q[0];
rz(2.1513596) q[1];
sx q[1];
rz(-2.4001887) q[1];
sx q[1];
rz(2.5726817) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.617482) q[0];
sx q[0];
rz(-2.4883399) q[0];
sx q[0];
rz(2.0481443) q[0];
x q[1];
rz(-2.8773035) q[2];
sx q[2];
rz(-1.2131872) q[2];
sx q[2];
rz(2.8823095) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68142707) q[1];
sx q[1];
rz(-0.54124628) q[1];
sx q[1];
rz(0.98249225) q[1];
rz(2.2248771) q[3];
sx q[3];
rz(-2.205483) q[3];
sx q[3];
rz(-3.0460639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27900055) q[2];
sx q[2];
rz(-2.180884) q[2];
sx q[2];
rz(-0.29435364) q[2];
rz(0.66458464) q[3];
sx q[3];
rz(-2.3817101) q[3];
sx q[3];
rz(-1.4544646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45254529) q[0];
sx q[0];
rz(-1.4280467) q[0];
sx q[0];
rz(2.8277165) q[0];
rz(2.2398056) q[1];
sx q[1];
rz(-1.7867463) q[1];
sx q[1];
rz(-1.4656167) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86014245) q[0];
sx q[0];
rz(-1.4518132) q[0];
sx q[0];
rz(-3.1156179) q[0];
x q[1];
rz(1.9889262) q[2];
sx q[2];
rz(-0.2731495) q[2];
sx q[2];
rz(1.7371617) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3180518) q[1];
sx q[1];
rz(-2.7465804) q[1];
sx q[1];
rz(-2.945963) q[1];
rz(-1.5336897) q[3];
sx q[3];
rz(-0.94054619) q[3];
sx q[3];
rz(-2.700151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19447154) q[2];
sx q[2];
rz(-0.7496382) q[2];
sx q[2];
rz(-2.0571902) q[2];
rz(-0.32367745) q[3];
sx q[3];
rz(-0.71666986) q[3];
sx q[3];
rz(-0.22423854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9913469) q[0];
sx q[0];
rz(-2.8960189) q[0];
sx q[0];
rz(-2.7299951) q[0];
rz(-2.4161074) q[1];
sx q[1];
rz(-1.2440224) q[1];
sx q[1];
rz(-1.7405608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66564225) q[0];
sx q[0];
rz(-2.4635297) q[0];
sx q[0];
rz(0.31291385) q[0];
x q[1];
rz(-0.13757041) q[2];
sx q[2];
rz(-2.3964632) q[2];
sx q[2];
rz(-2.957475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8795149) q[1];
sx q[1];
rz(-1.8873155) q[1];
sx q[1];
rz(-2.7807117) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15604354) q[3];
sx q[3];
rz(-0.45609176) q[3];
sx q[3];
rz(1.9182916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.97658849) q[2];
sx q[2];
rz(-0.62825957) q[2];
sx q[2];
rz(-1.7151625) q[2];
rz(0.12380883) q[3];
sx q[3];
rz(-1.0694458) q[3];
sx q[3];
rz(-0.4933221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2774778) q[0];
sx q[0];
rz(-1.9954229) q[0];
sx q[0];
rz(-1.8051099) q[0];
rz(-2.2026964) q[1];
sx q[1];
rz(-1.122033) q[1];
sx q[1];
rz(-0.28405651) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1079252) q[0];
sx q[0];
rz(-1.2311651) q[0];
sx q[0];
rz(-2.4976298) q[0];
rz(-pi) q[1];
rz(-2.192382) q[2];
sx q[2];
rz(-0.27602592) q[2];
sx q[2];
rz(-0.29601184) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1394609) q[1];
sx q[1];
rz(-1.5393917) q[1];
sx q[1];
rz(-3.1314) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18064336) q[3];
sx q[3];
rz(-1.7076645) q[3];
sx q[3];
rz(-1.6302072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17860086) q[2];
sx q[2];
rz(-2.4602349) q[2];
sx q[2];
rz(-2.2006688) q[2];
rz(-3.1407147) q[3];
sx q[3];
rz(-1.2255171) q[3];
sx q[3];
rz(3.0555449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0451999) q[0];
sx q[0];
rz(-0.89458507) q[0];
sx q[0];
rz(-2.9659502) q[0];
rz(1.9215709) q[1];
sx q[1];
rz(-2.2717387) q[1];
sx q[1];
rz(-2.2697935) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5822495) q[0];
sx q[0];
rz(-1.7212369) q[0];
sx q[0];
rz(-1.3036672) q[0];
rz(-pi) q[1];
rz(0.37105889) q[2];
sx q[2];
rz(-2.3239229) q[2];
sx q[2];
rz(1.0196067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2928472) q[1];
sx q[1];
rz(-2.1715144) q[1];
sx q[1];
rz(1.377618) q[1];
rz(-pi) q[2];
rz(-0.11856793) q[3];
sx q[3];
rz(-2.5817079) q[3];
sx q[3];
rz(0.15794755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.72622) q[2];
sx q[2];
rz(-1.8612334) q[2];
sx q[2];
rz(1.0225164) q[2];
rz(0.38698777) q[3];
sx q[3];
rz(-1.756668) q[3];
sx q[3];
rz(-0.81138396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33685327) q[0];
sx q[0];
rz(-1.9751208) q[0];
sx q[0];
rz(-2.8785896) q[0];
rz(-0.31245843) q[1];
sx q[1];
rz(-2.0461021) q[1];
sx q[1];
rz(1.9591263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4561037) q[0];
sx q[0];
rz(-2.0747599) q[0];
sx q[0];
rz(-2.9188372) q[0];
rz(-pi) q[1];
rz(-0.83177213) q[2];
sx q[2];
rz(-1.1397994) q[2];
sx q[2];
rz(0.43048358) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2559465) q[1];
sx q[1];
rz(-1.3722536) q[1];
sx q[1];
rz(-2.7629653) q[1];
rz(1.0013323) q[3];
sx q[3];
rz(-2.0234442) q[3];
sx q[3];
rz(2.4695726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.075228127) q[2];
sx q[2];
rz(-1.790739) q[2];
sx q[2];
rz(2.8933375) q[2];
rz(2.9634641) q[3];
sx q[3];
rz(-2.0592561) q[3];
sx q[3];
rz(-2.3586912) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7149776) q[0];
sx q[0];
rz(-2.7474032) q[0];
sx q[0];
rz(2.07975) q[0];
rz(1.3007523) q[1];
sx q[1];
rz(-1.5251093) q[1];
sx q[1];
rz(1.1311857) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12515629) q[0];
sx q[0];
rz(-1.1068692) q[0];
sx q[0];
rz(2.6317276) q[0];
x q[1];
rz(0.71604095) q[2];
sx q[2];
rz(-0.63010213) q[2];
sx q[2];
rz(2.1269456) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2857745) q[1];
sx q[1];
rz(-2.7322391) q[1];
sx q[1];
rz(-0.23366433) q[1];
rz(-pi) q[2];
rz(0.33343306) q[3];
sx q[3];
rz(-2.5134183) q[3];
sx q[3];
rz(-1.4905358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0807557) q[2];
sx q[2];
rz(-2.8751825) q[2];
sx q[2];
rz(0.60144919) q[2];
rz(-2.1001749) q[3];
sx q[3];
rz(-1.4789378) q[3];
sx q[3];
rz(-1.7634332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1489442) q[0];
sx q[0];
rz(-2.5288378) q[0];
sx q[0];
rz(0.80768325) q[0];
rz(3.0638937) q[1];
sx q[1];
rz(-2.6014889) q[1];
sx q[1];
rz(-1.4059975) q[1];
rz(2.7355657) q[2];
sx q[2];
rz(-2.1320504) q[2];
sx q[2];
rz(-0.8036094) q[2];
rz(0.64114707) q[3];
sx q[3];
rz(-2.2401886) q[3];
sx q[3];
rz(-0.2284579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
