OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90836877) q[0];
sx q[0];
rz(4.6908661) q[0];
sx q[0];
rz(9.6022845) q[0];
rz(-1.5818051) q[1];
sx q[1];
rz(-1.3270562) q[1];
sx q[1];
rz(-1.9628734) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66087276) q[0];
sx q[0];
rz(-2.5205223) q[0];
sx q[0];
rz(1.2476515) q[0];
rz(-pi) q[1];
rz(-0.92843797) q[2];
sx q[2];
rz(-1.3283159) q[2];
sx q[2];
rz(1.2280304) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.56621579) q[1];
sx q[1];
rz(-1.0027835) q[1];
sx q[1];
rz(-2.6188208) q[1];
rz(1.7186382) q[3];
sx q[3];
rz(-2.4018722) q[3];
sx q[3];
rz(-1.5613558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.061964758) q[2];
sx q[2];
rz(-0.92755932) q[2];
sx q[2];
rz(2.8947042) q[2];
rz(-2.7039792) q[3];
sx q[3];
rz(-2.103431) q[3];
sx q[3];
rz(2.8810697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058411773) q[0];
sx q[0];
rz(-2.658168) q[0];
sx q[0];
rz(-2.0368982) q[0];
rz(-0.80766922) q[1];
sx q[1];
rz(-0.75391114) q[1];
sx q[1];
rz(-2.5025936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4527567) q[0];
sx q[0];
rz(-0.3301183) q[0];
sx q[0];
rz(0.56906505) q[0];
rz(-2.866687) q[2];
sx q[2];
rz(-1.5228809) q[2];
sx q[2];
rz(-0.78457441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2874291) q[1];
sx q[1];
rz(-0.87814444) q[1];
sx q[1];
rz(-2.0504897) q[1];
x q[2];
rz(0.2868594) q[3];
sx q[3];
rz(-1.8897227) q[3];
sx q[3];
rz(0.87251679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.69940007) q[2];
sx q[2];
rz(-1.7510119) q[2];
sx q[2];
rz(0.81713444) q[2];
rz(-0.52554956) q[3];
sx q[3];
rz(-2.3216129) q[3];
sx q[3];
rz(1.9513963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84080559) q[0];
sx q[0];
rz(-1.5910933) q[0];
sx q[0];
rz(2.9817885) q[0];
rz(-2.6218759) q[1];
sx q[1];
rz(-2.0964041) q[1];
sx q[1];
rz(2.2249178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0063112886) q[0];
sx q[0];
rz(-2.1410628) q[0];
sx q[0];
rz(1.6018014) q[0];
x q[1];
rz(-1.3266852) q[2];
sx q[2];
rz(-2.296431) q[2];
sx q[2];
rz(-0.47032088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8868272) q[1];
sx q[1];
rz(-1.5228094) q[1];
sx q[1];
rz(1.3814203) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4363719) q[3];
sx q[3];
rz(-2.1791576) q[3];
sx q[3];
rz(0.20316635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2436169) q[2];
sx q[2];
rz(-0.50668442) q[2];
sx q[2];
rz(1.766073) q[2];
rz(-3.1084133) q[3];
sx q[3];
rz(-1.5539919) q[3];
sx q[3];
rz(-0.83056617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631977) q[0];
sx q[0];
rz(-2.4802408) q[0];
sx q[0];
rz(2.1589494) q[0];
rz(-2.4962418) q[1];
sx q[1];
rz(-1.2256365) q[1];
sx q[1];
rz(2.7243848) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.062995) q[0];
sx q[0];
rz(-1.0630442) q[0];
sx q[0];
rz(1.9168454) q[0];
rz(-1.7424612) q[2];
sx q[2];
rz(-2.7560803) q[2];
sx q[2];
rz(-2.3814892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9703435) q[1];
sx q[1];
rz(-1.6303195) q[1];
sx q[1];
rz(-0.67386676) q[1];
rz(-0.44261114) q[3];
sx q[3];
rz(-2.5472982) q[3];
sx q[3];
rz(1.8123008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8365606) q[2];
sx q[2];
rz(-2.8053668) q[2];
sx q[2];
rz(-1.9545371) q[2];
rz(0.120397) q[3];
sx q[3];
rz(-1.7359366) q[3];
sx q[3];
rz(0.19336893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8315941) q[0];
sx q[0];
rz(-3.031142) q[0];
sx q[0];
rz(2.4647392) q[0];
rz(1.6090769) q[1];
sx q[1];
rz(-1.1108421) q[1];
sx q[1];
rz(0.19930509) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4368962) q[0];
sx q[0];
rz(-1.4802209) q[0];
sx q[0];
rz(-1.8245839) q[0];
rz(-1.4485613) q[2];
sx q[2];
rz(-1.9150371) q[2];
sx q[2];
rz(-0.93512541) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2414294) q[1];
sx q[1];
rz(-1.8668264) q[1];
sx q[1];
rz(2.7524292) q[1];
x q[2];
rz(1.293888) q[3];
sx q[3];
rz(-1.5770638) q[3];
sx q[3];
rz(2.4771877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4946332) q[2];
sx q[2];
rz(-0.55611742) q[2];
sx q[2];
rz(0.49590084) q[2];
rz(2.2206709) q[3];
sx q[3];
rz(-2.0995188) q[3];
sx q[3];
rz(-2.8310149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.6344236) q[0];
sx q[0];
rz(-1.2260219) q[0];
sx q[0];
rz(-0.41967151) q[0];
rz(-0.93027973) q[1];
sx q[1];
rz(-1.0134965) q[1];
sx q[1];
rz(1.0118265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7714846) q[0];
sx q[0];
rz(-0.91910942) q[0];
sx q[0];
rz(-1.6328567) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65134766) q[2];
sx q[2];
rz(-0.22964165) q[2];
sx q[2];
rz(2.1852213) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7020445) q[1];
sx q[1];
rz(-1.7782033) q[1];
sx q[1];
rz(-1.009936) q[1];
x q[2];
rz(-2.9568372) q[3];
sx q[3];
rz(-2.1596381) q[3];
sx q[3];
rz(-2.7240804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73152995) q[2];
sx q[2];
rz(-0.56649929) q[2];
sx q[2];
rz(-2.0274963) q[2];
rz(1.7690432) q[3];
sx q[3];
rz(-1.0053582) q[3];
sx q[3];
rz(-2.4323288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0396742) q[0];
sx q[0];
rz(-1.2693951) q[0];
sx q[0];
rz(-1.4920379) q[0];
rz(-2.6986625) q[1];
sx q[1];
rz(-1.2716753) q[1];
sx q[1];
rz(2.2992004) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32552606) q[0];
sx q[0];
rz(-0.93579817) q[0];
sx q[0];
rz(0.86845894) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2350253) q[2];
sx q[2];
rz(-1.3593055) q[2];
sx q[2];
rz(-1.4480526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2017224) q[1];
sx q[1];
rz(-2.2186167) q[1];
sx q[1];
rz(2.617791) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.962389) q[3];
sx q[3];
rz(-0.11193724) q[3];
sx q[3];
rz(-0.051210784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3288021) q[2];
sx q[2];
rz(-0.91965681) q[2];
sx q[2];
rz(0.021050464) q[2];
rz(-2.1444495) q[3];
sx q[3];
rz(-0.54233426) q[3];
sx q[3];
rz(1.5525345) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7954623) q[0];
sx q[0];
rz(-1.2393476) q[0];
sx q[0];
rz(-1.963266) q[0];
rz(1.0109674) q[1];
sx q[1];
rz(-1.4821472) q[1];
sx q[1];
rz(-0.3709901) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8218433) q[0];
sx q[0];
rz(-0.25938636) q[0];
sx q[0];
rz(-0.46617561) q[0];
rz(-pi) q[1];
rz(-1.0576789) q[2];
sx q[2];
rz(-2.9484595) q[2];
sx q[2];
rz(0.73207515) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4996262) q[1];
sx q[1];
rz(-1.8607983) q[1];
sx q[1];
rz(-0.15695842) q[1];
x q[2];
rz(2.8274687) q[3];
sx q[3];
rz(-1.1289136) q[3];
sx q[3];
rz(-1.8754547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3498938) q[2];
sx q[2];
rz(-1.0827622) q[2];
sx q[2];
rz(0.91378158) q[2];
rz(-2.6036116) q[3];
sx q[3];
rz(-2.3086083) q[3];
sx q[3];
rz(0.37916455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0898042) q[0];
sx q[0];
rz(-1.3642949) q[0];
sx q[0];
rz(-0.23189932) q[0];
rz(-1.505835) q[1];
sx q[1];
rz(-1.4488522) q[1];
sx q[1];
rz(-2.4969782) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1132677) q[0];
sx q[0];
rz(-1.2282983) q[0];
sx q[0];
rz(1.2069846) q[0];
rz(-3.1050131) q[2];
sx q[2];
rz(-2.5365005) q[2];
sx q[2];
rz(-0.99672752) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48460173) q[1];
sx q[1];
rz(-1.4047523) q[1];
sx q[1];
rz(0.019980494) q[1];
rz(-1.1357952) q[3];
sx q[3];
rz(-0.91640546) q[3];
sx q[3];
rz(-1.1799348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5171648) q[2];
sx q[2];
rz(-1.5394053) q[2];
sx q[2];
rz(1.2702764) q[2];
rz(-2.7020448) q[3];
sx q[3];
rz(-2.5054273) q[3];
sx q[3];
rz(2.6308681) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9257833) q[0];
sx q[0];
rz(-1.2405115) q[0];
sx q[0];
rz(-2.1155758) q[0];
rz(0.45694524) q[1];
sx q[1];
rz(-0.89070717) q[1];
sx q[1];
rz(2.2073726) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81722014) q[0];
sx q[0];
rz(-1.1486832) q[0];
sx q[0];
rz(-1.425231) q[0];
rz(0.1100357) q[2];
sx q[2];
rz(-1.8241183) q[2];
sx q[2];
rz(0.91722721) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9125713) q[1];
sx q[1];
rz(-2.3816436) q[1];
sx q[1];
rz(0.2473803) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85984965) q[3];
sx q[3];
rz(-0.78892869) q[3];
sx q[3];
rz(-3.0956059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43941867) q[2];
sx q[2];
rz(-2.4472523) q[2];
sx q[2];
rz(1.0609421) q[2];
rz(0.98961467) q[3];
sx q[3];
rz(-1.7695534) q[3];
sx q[3];
rz(2.8770679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92985741) q[0];
sx q[0];
rz(-2.0462357) q[0];
sx q[0];
rz(-0.76475058) q[0];
rz(1.1872956) q[1];
sx q[1];
rz(-1.7056414) q[1];
sx q[1];
rz(2.5797226) q[1];
rz(-1.622772) q[2];
sx q[2];
rz(-1.5511654) q[2];
sx q[2];
rz(1.8775107) q[2];
rz(0.46502385) q[3];
sx q[3];
rz(-1.1661543) q[3];
sx q[3];
rz(2.3247624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
