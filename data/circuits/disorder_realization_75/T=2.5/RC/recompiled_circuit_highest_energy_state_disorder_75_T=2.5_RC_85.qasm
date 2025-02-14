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
rz(1.6504352) q[0];
sx q[0];
rz(-0.37547922) q[0];
sx q[0];
rz(-0.39146358) q[0];
rz(-2.904881) q[1];
sx q[1];
rz(-1.038329) q[1];
sx q[1];
rz(2.2957323) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5206155) q[0];
sx q[0];
rz(-1.3987887) q[0];
sx q[0];
rz(-0.67396236) q[0];
rz(-pi) q[1];
rz(1.1089788) q[2];
sx q[2];
rz(-1.9792288) q[2];
sx q[2];
rz(1.3549182) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9698025) q[1];
sx q[1];
rz(-0.39229052) q[1];
sx q[1];
rz(-0.7902625) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0580642) q[3];
sx q[3];
rz(-1.8780276) q[3];
sx q[3];
rz(1.787825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.55964959) q[2];
sx q[2];
rz(-2.1008284) q[2];
sx q[2];
rz(-2.9259658) q[2];
rz(-0.98254472) q[3];
sx q[3];
rz(-1.6590786) q[3];
sx q[3];
rz(1.8854347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243645) q[0];
sx q[0];
rz(-0.094014458) q[0];
sx q[0];
rz(2.7650058) q[0];
rz(1.6603893) q[1];
sx q[1];
rz(-1.1777638) q[1];
sx q[1];
rz(1.0053763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1373077) q[0];
sx q[0];
rz(-1.8958603) q[0];
sx q[0];
rz(0.0079197366) q[0];
rz(-pi) q[1];
rz(2.2475904) q[2];
sx q[2];
rz(-1.4419705) q[2];
sx q[2];
rz(-1.027838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.95208827) q[1];
sx q[1];
rz(-0.90248195) q[1];
sx q[1];
rz(-1.52723) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5672417) q[3];
sx q[3];
rz(-0.6308517) q[3];
sx q[3];
rz(1.3441946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5542095) q[2];
sx q[2];
rz(-0.46450928) q[2];
sx q[2];
rz(-1.0760388) q[2];
rz(0.46766034) q[3];
sx q[3];
rz(-2.3105919) q[3];
sx q[3];
rz(-0.20461288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6429546) q[0];
sx q[0];
rz(-2.0641646) q[0];
sx q[0];
rz(1.1109918) q[0];
rz(-2.8338762) q[1];
sx q[1];
rz(-1.6650763) q[1];
sx q[1];
rz(-1.841338) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4487023) q[0];
sx q[0];
rz(-1.4822791) q[0];
sx q[0];
rz(0.55175106) q[0];
x q[1];
rz(0.05540313) q[2];
sx q[2];
rz(-0.57864648) q[2];
sx q[2];
rz(2.5134653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.17323) q[1];
sx q[1];
rz(-0.97665962) q[1];
sx q[1];
rz(2.6244375) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6407317) q[3];
sx q[3];
rz(-1.259049) q[3];
sx q[3];
rz(0.65697981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51851455) q[2];
sx q[2];
rz(-1.0479505) q[2];
sx q[2];
rz(2.9497362) q[2];
rz(2.0375371) q[3];
sx q[3];
rz(-0.42686978) q[3];
sx q[3];
rz(-1.2087315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.0693102) q[0];
sx q[0];
rz(-0.59232124) q[0];
sx q[0];
rz(0.88822547) q[0];
rz(0.2725254) q[1];
sx q[1];
rz(-1.2963632) q[1];
sx q[1];
rz(1.9857508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4711535) q[0];
sx q[0];
rz(-1.5799755) q[0];
sx q[0];
rz(1.2769475) q[0];
rz(-1.5290909) q[2];
sx q[2];
rz(-1.0805849) q[2];
sx q[2];
rz(0.34402564) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9458876) q[1];
sx q[1];
rz(-0.81521704) q[1];
sx q[1];
rz(-2.7570711) q[1];
x q[2];
rz(-0.0083752092) q[3];
sx q[3];
rz(-2.6937458) q[3];
sx q[3];
rz(0.26219246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83879519) q[2];
sx q[2];
rz(-0.28442997) q[2];
sx q[2];
rz(-2.3885041) q[2];
rz(-0.49790844) q[3];
sx q[3];
rz(-1.4830736) q[3];
sx q[3];
rz(-0.30043093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.65557757) q[0];
sx q[0];
rz(-0.97633728) q[0];
sx q[0];
rz(1.2855541) q[0];
rz(-1.1898419) q[1];
sx q[1];
rz(-0.68584502) q[1];
sx q[1];
rz(1.5815585) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043997633) q[0];
sx q[0];
rz(-1.4420995) q[0];
sx q[0];
rz(3.0630847) q[0];
rz(-pi) q[1];
rz(2.8578256) q[2];
sx q[2];
rz(-1.7161676) q[2];
sx q[2];
rz(-1.870537) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6717357) q[1];
sx q[1];
rz(-1.3570045) q[1];
sx q[1];
rz(0.44624568) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64636773) q[3];
sx q[3];
rz(-2.8146675) q[3];
sx q[3];
rz(-1.2267139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.84216422) q[2];
sx q[2];
rz(-2.2101768) q[2];
sx q[2];
rz(0.13775873) q[2];
rz(2.6970741) q[3];
sx q[3];
rz(-1.7701745) q[3];
sx q[3];
rz(0.81454149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21823068) q[0];
sx q[0];
rz(-0.87319279) q[0];
sx q[0];
rz(-2.5905304) q[0];
rz(-0.29523826) q[1];
sx q[1];
rz(-2.4199838) q[1];
sx q[1];
rz(2.3982184) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0449042) q[0];
sx q[0];
rz(-1.4375234) q[0];
sx q[0];
rz(-3.0267984) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34140519) q[2];
sx q[2];
rz(-1.7936754) q[2];
sx q[2];
rz(2.6472732) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.36657143) q[1];
sx q[1];
rz(-1.4719392) q[1];
sx q[1];
rz(-2.2926383) q[1];
x q[2];
rz(2.1990699) q[3];
sx q[3];
rz(-0.97369558) q[3];
sx q[3];
rz(-1.7973078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4117406) q[2];
sx q[2];
rz(-2.9509632) q[2];
sx q[2];
rz(-0.28212485) q[2];
rz(0.43635803) q[3];
sx q[3];
rz(-1.8264419) q[3];
sx q[3];
rz(-0.25562975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.2082763) q[0];
sx q[0];
rz(-2.0774807) q[0];
sx q[0];
rz(-2.5489885) q[0];
rz(-1.1366049) q[1];
sx q[1];
rz(-1.2717783) q[1];
sx q[1];
rz(2.073854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.655433) q[0];
sx q[0];
rz(-0.91595931) q[0];
sx q[0];
rz(-1.7688807) q[0];
x q[1];
rz(1.4456533) q[2];
sx q[2];
rz(-1.2045025) q[2];
sx q[2];
rz(2.3082781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6954036) q[1];
sx q[1];
rz(-1.0985038) q[1];
sx q[1];
rz(2.8279734) q[1];
x q[2];
rz(-0.41530825) q[3];
sx q[3];
rz(-1.7303932) q[3];
sx q[3];
rz(0.94091614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9390255) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(-0.97958809) q[2];
rz(0.12796399) q[3];
sx q[3];
rz(-2.1991859) q[3];
sx q[3];
rz(1.452182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99280438) q[0];
sx q[0];
rz(-2.1908741) q[0];
sx q[0];
rz(3.0606781) q[0];
rz(3.0203536) q[1];
sx q[1];
rz(-0.67758766) q[1];
sx q[1];
rz(-1.1239207) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2291717) q[0];
sx q[0];
rz(-1.9397368) q[0];
sx q[0];
rz(-2.3652618) q[0];
rz(-pi) q[1];
rz(-3.0380855) q[2];
sx q[2];
rz(-0.40914224) q[2];
sx q[2];
rz(1.4827267) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2595265) q[1];
sx q[1];
rz(-1.2203382) q[1];
sx q[1];
rz(1.4653852) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0392022) q[3];
sx q[3];
rz(-1.7178159) q[3];
sx q[3];
rz(-1.3519629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.27440444) q[2];
sx q[2];
rz(-0.97325456) q[2];
sx q[2];
rz(-1.6483866) q[2];
rz(0.46227208) q[3];
sx q[3];
rz(-0.8148163) q[3];
sx q[3];
rz(1.7260684) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129775) q[0];
sx q[0];
rz(-1.3884437) q[0];
sx q[0];
rz(-1.9762565) q[0];
rz(-0.041821592) q[1];
sx q[1];
rz(-0.98894293) q[1];
sx q[1];
rz(-1.3335386) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8671252) q[0];
sx q[0];
rz(-2.2750912) q[0];
sx q[0];
rz(2.254266) q[0];
rz(1.5666577) q[2];
sx q[2];
rz(-1.2519826) q[2];
sx q[2];
rz(-0.37296527) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11631498) q[1];
sx q[1];
rz(-1.9851369) q[1];
sx q[1];
rz(-1.6700909) q[1];
x q[2];
rz(1.4708552) q[3];
sx q[3];
rz(-1.458711) q[3];
sx q[3];
rz(-2.374927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24451438) q[2];
sx q[2];
rz(-1.3099193) q[2];
sx q[2];
rz(2.6778527) q[2];
rz(2.5044299) q[3];
sx q[3];
rz(-1.623268) q[3];
sx q[3];
rz(-2.7630828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31507444) q[0];
sx q[0];
rz(-1.0239064) q[0];
sx q[0];
rz(1.4950604) q[0];
rz(0.82978326) q[1];
sx q[1];
rz(-1.2636431) q[1];
sx q[1];
rz(-1.319938) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76413122) q[0];
sx q[0];
rz(-0.39239663) q[0];
sx q[0];
rz(-2.8026587) q[0];
rz(1.3426724) q[2];
sx q[2];
rz(-1.14883) q[2];
sx q[2];
rz(-2.162148) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1935008) q[1];
sx q[1];
rz(-1.337916) q[1];
sx q[1];
rz(1.3629254) q[1];
rz(2.8828524) q[3];
sx q[3];
rz(-1.0240533) q[3];
sx q[3];
rz(-0.26165831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39184555) q[2];
sx q[2];
rz(-1.7228935) q[2];
sx q[2];
rz(0.5113655) q[2];
rz(0.5558719) q[3];
sx q[3];
rz(-2.0248196) q[3];
sx q[3];
rz(-2.6218124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2340672) q[0];
sx q[0];
rz(-3.0341442) q[0];
sx q[0];
rz(-2.7035614) q[0];
rz(0.058649339) q[1];
sx q[1];
rz(-1.5025243) q[1];
sx q[1];
rz(2.0741838) q[1];
rz(1.5127586) q[2];
sx q[2];
rz(-2.2709202) q[2];
sx q[2];
rz(1.5135598) q[2];
rz(-2.5142148) q[3];
sx q[3];
rz(-1.3199721) q[3];
sx q[3];
rz(-1.6047431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
