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
rz(-3.0523678) q[0];
sx q[0];
rz(-2.3951946) q[0];
sx q[0];
rz(2.449692) q[0];
rz(0.46299419) q[1];
sx q[1];
rz(-1.6166592) q[1];
sx q[1];
rz(-1.0301627) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2041516) q[0];
sx q[0];
rz(-0.61249706) q[0];
sx q[0];
rz(-0.4093616) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9691508) q[2];
sx q[2];
rz(-1.4015369) q[2];
sx q[2];
rz(-1.5072849) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4575689) q[1];
sx q[1];
rz(-1.6872703) q[1];
sx q[1];
rz(1.2435786) q[1];
rz(-pi) q[2];
x q[2];
rz(2.452677) q[3];
sx q[3];
rz(-1.6647881) q[3];
sx q[3];
rz(1.1979562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2618711) q[2];
sx q[2];
rz(-2.3679831) q[2];
sx q[2];
rz(0.36641463) q[2];
rz(-2.053818) q[3];
sx q[3];
rz(-2.352114) q[3];
sx q[3];
rz(-2.4521949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.825603) q[0];
sx q[0];
rz(-3.0593384) q[0];
sx q[0];
rz(-0.91973037) q[0];
rz(-0.74291825) q[1];
sx q[1];
rz(-1.2580322) q[1];
sx q[1];
rz(2.0222576) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9237083) q[0];
sx q[0];
rz(-0.15734921) q[0];
sx q[0];
rz(1.1697392) q[0];
x q[1];
rz(-0.2074457) q[2];
sx q[2];
rz(-2.4118057) q[2];
sx q[2];
rz(1.6285673) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7537958) q[1];
sx q[1];
rz(-2.2860074) q[1];
sx q[1];
rz(-1.254735) q[1];
x q[2];
rz(-2.2833586) q[3];
sx q[3];
rz(-2.4284114) q[3];
sx q[3];
rz(-0.76480097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4982345) q[2];
sx q[2];
rz(-1.9469807) q[2];
sx q[2];
rz(-3.0894792) q[2];
rz(1.107996) q[3];
sx q[3];
rz(-2.2823157) q[3];
sx q[3];
rz(3.0767379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4614918) q[0];
sx q[0];
rz(-1.2315741) q[0];
sx q[0];
rz(-1.8517866) q[0];
rz(-2.3884804) q[1];
sx q[1];
rz(-1.6878004) q[1];
sx q[1];
rz(2.6444816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4418934) q[0];
sx q[0];
rz(-1.178368) q[0];
sx q[0];
rz(2.5470069) q[0];
rz(2.3894823) q[2];
sx q[2];
rz(-2.3083938) q[2];
sx q[2];
rz(-1.1583661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.048084413) q[1];
sx q[1];
rz(-1.1968379) q[1];
sx q[1];
rz(-1.5300586) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6670973) q[3];
sx q[3];
rz(-2.2212006) q[3];
sx q[3];
rz(-1.4092829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6734753) q[2];
sx q[2];
rz(-0.27286369) q[2];
sx q[2];
rz(2.4172778) q[2];
rz(-2.916548) q[3];
sx q[3];
rz(-1.2757755) q[3];
sx q[3];
rz(0.79506522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.986079) q[0];
sx q[0];
rz(-2.2561095) q[0];
sx q[0];
rz(0.28050637) q[0];
rz(1.4836503) q[1];
sx q[1];
rz(-1.2018485) q[1];
sx q[1];
rz(2.6703506) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12898937) q[0];
sx q[0];
rz(-1.5480729) q[0];
sx q[0];
rz(1.9381592) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.066023237) q[2];
sx q[2];
rz(-1.4581265) q[2];
sx q[2];
rz(-0.45399912) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.60741456) q[1];
sx q[1];
rz(-1.3001658) q[1];
sx q[1];
rz(-3.1128008) q[1];
rz(0.56994168) q[3];
sx q[3];
rz(-2.2472858) q[3];
sx q[3];
rz(-2.3968056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9435141) q[2];
sx q[2];
rz(-2.5835865) q[2];
sx q[2];
rz(1.0151218) q[2];
rz(0.73457581) q[3];
sx q[3];
rz(-1.3646804) q[3];
sx q[3];
rz(-2.6034897) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15767844) q[0];
sx q[0];
rz(-1.2175125) q[0];
sx q[0];
rz(2.1814046) q[0];
rz(0.045711191) q[1];
sx q[1];
rz(-0.8784596) q[1];
sx q[1];
rz(-0.033230573) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93532256) q[0];
sx q[0];
rz(-2.0747796) q[0];
sx q[0];
rz(-2.5978243) q[0];
rz(1.5157891) q[2];
sx q[2];
rz(-0.26794456) q[2];
sx q[2];
rz(-2.1017113) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3551801) q[1];
sx q[1];
rz(-2.7828453) q[1];
sx q[1];
rz(-1.1915818) q[1];
x q[2];
rz(1.9779677) q[3];
sx q[3];
rz(-1.1615175) q[3];
sx q[3];
rz(0.57151505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9783322) q[2];
sx q[2];
rz(-2.1425118) q[2];
sx q[2];
rz(-0.31326374) q[2];
rz(2.7836109) q[3];
sx q[3];
rz(-2.101818) q[3];
sx q[3];
rz(0.09859214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1545496) q[0];
sx q[0];
rz(-1.6313666) q[0];
sx q[0];
rz(-0.48598591) q[0];
rz(-1.9588574) q[1];
sx q[1];
rz(-0.69570884) q[1];
sx q[1];
rz(0.33014306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3957452) q[0];
sx q[0];
rz(-1.2735277) q[0];
sx q[0];
rz(-0.031121522) q[0];
rz(-pi) q[1];
rz(0.62689828) q[2];
sx q[2];
rz(-0.90993249) q[2];
sx q[2];
rz(-0.063274296) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9818519) q[1];
sx q[1];
rz(-1.0969437) q[1];
sx q[1];
rz(-0.57094806) q[1];
rz(2.7948912) q[3];
sx q[3];
rz(-1.3606405) q[3];
sx q[3];
rz(0.29276785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5064454) q[2];
sx q[2];
rz(-2.5250285) q[2];
sx q[2];
rz(-2.7849624) q[2];
rz(0.59456524) q[3];
sx q[3];
rz(-1.6584572) q[3];
sx q[3];
rz(2.4026292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67721382) q[0];
sx q[0];
rz(-0.76736275) q[0];
sx q[0];
rz(-2.5546524) q[0];
rz(0.37239536) q[1];
sx q[1];
rz(-1.3601235) q[1];
sx q[1];
rz(2.8940103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0599759) q[0];
sx q[0];
rz(-0.46591407) q[0];
sx q[0];
rz(1.5008677) q[0];
x q[1];
rz(1.814194) q[2];
sx q[2];
rz(-2.1888424) q[2];
sx q[2];
rz(-0.40153654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0928318) q[1];
sx q[1];
rz(-2.2032461) q[1];
sx q[1];
rz(-0.41776843) q[1];
x q[2];
rz(1.8195175) q[3];
sx q[3];
rz(-1.910672) q[3];
sx q[3];
rz(-2.9074977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8783012) q[2];
sx q[2];
rz(-1.2211439) q[2];
sx q[2];
rz(-2.9325874) q[2];
rz(2.5448997) q[3];
sx q[3];
rz(-2.7946819) q[3];
sx q[3];
rz(-1.5359115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40808943) q[0];
sx q[0];
rz(-0.96036378) q[0];
sx q[0];
rz(-2.7485513) q[0];
rz(2.4709002) q[1];
sx q[1];
rz(-1.9978943) q[1];
sx q[1];
rz(-1.447698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2358952) q[0];
sx q[0];
rz(-1.4465155) q[0];
sx q[0];
rz(-2.0688562) q[0];
rz(-pi) q[1];
rz(-0.31140458) q[2];
sx q[2];
rz(-0.97955737) q[2];
sx q[2];
rz(-0.55252749) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5504584) q[1];
sx q[1];
rz(-1.0626284) q[1];
sx q[1];
rz(-0.9891467) q[1];
rz(-pi) q[2];
rz(2.0441229) q[3];
sx q[3];
rz(-2.3103263) q[3];
sx q[3];
rz(-2.748718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4286917) q[2];
sx q[2];
rz(-2.28076) q[2];
sx q[2];
rz(-0.11670308) q[2];
rz(2.1278925) q[3];
sx q[3];
rz(-1.2319535) q[3];
sx q[3];
rz(-1.8359756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75334615) q[0];
sx q[0];
rz(-1.6755063) q[0];
sx q[0];
rz(-1.7062794) q[0];
rz(-2.1460136) q[1];
sx q[1];
rz(-1.8228056) q[1];
sx q[1];
rz(-2.6011655) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0198959) q[0];
sx q[0];
rz(-1.3084687) q[0];
sx q[0];
rz(0.4957677) q[0];
x q[1];
rz(-2.5535271) q[2];
sx q[2];
rz(-0.99074513) q[2];
sx q[2];
rz(1.560245) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82481532) q[1];
sx q[1];
rz(-1.1939923) q[1];
sx q[1];
rz(2.3945859) q[1];
rz(-1.1630658) q[3];
sx q[3];
rz(-0.98568788) q[3];
sx q[3];
rz(-1.9118904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.871375) q[2];
sx q[2];
rz(-1.3048708) q[2];
sx q[2];
rz(0.77667856) q[2];
rz(-0.03260472) q[3];
sx q[3];
rz(-1.6083345) q[3];
sx q[3];
rz(-1.5743272) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33083522) q[0];
sx q[0];
rz(-0.36892712) q[0];
sx q[0];
rz(-0.43357968) q[0];
rz(-2.4724204) q[1];
sx q[1];
rz(-1.9586261) q[1];
sx q[1];
rz(2.7152667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0678372) q[0];
sx q[0];
rz(-1.3569419) q[0];
sx q[0];
rz(0.20238371) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2157006) q[2];
sx q[2];
rz(-0.85596687) q[2];
sx q[2];
rz(2.1260043) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59407104) q[1];
sx q[1];
rz(-1.921355) q[1];
sx q[1];
rz(-0.5188397) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4678585) q[3];
sx q[3];
rz(-0.70286432) q[3];
sx q[3];
rz(1.7029312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3907889) q[2];
sx q[2];
rz(-2.8533253) q[2];
sx q[2];
rz(-1.1061579) q[2];
rz(3.028051) q[3];
sx q[3];
rz(-2.2083486) q[3];
sx q[3];
rz(0.043702628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15585598) q[0];
sx q[0];
rz(-0.899213) q[0];
sx q[0];
rz(-0.56094299) q[0];
rz(-0.38324311) q[1];
sx q[1];
rz(-2.0189197) q[1];
sx q[1];
rz(-0.22519208) q[1];
rz(-2.6720895) q[2];
sx q[2];
rz(-0.60890763) q[2];
sx q[2];
rz(1.5694697) q[2];
rz(2.3229937) q[3];
sx q[3];
rz(-1.9055453) q[3];
sx q[3];
rz(0.26858134) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
