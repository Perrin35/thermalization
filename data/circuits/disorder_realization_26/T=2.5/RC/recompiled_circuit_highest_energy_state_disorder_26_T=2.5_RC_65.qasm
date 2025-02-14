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
rz(-2.2505724) q[0];
sx q[0];
rz(-1.855259) q[0];
sx q[0];
rz(2.2556055) q[0];
rz(-2.6180144) q[1];
sx q[1];
rz(-0.55966592) q[1];
sx q[1];
rz(1.3203415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65151764) q[0];
sx q[0];
rz(-2.1044255) q[0];
sx q[0];
rz(-0.20488157) q[0];
rz(-pi) q[1];
rz(0.34446005) q[2];
sx q[2];
rz(-2.3340324) q[2];
sx q[2];
rz(-1.5429614) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3229827) q[1];
sx q[1];
rz(-0.56639379) q[1];
sx q[1];
rz(-1.4569805) q[1];
rz(-pi) q[2];
rz(2.3778649) q[3];
sx q[3];
rz(-2.1565752) q[3];
sx q[3];
rz(1.8770742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3628799) q[2];
sx q[2];
rz(-1.3211297) q[2];
sx q[2];
rz(-2.2134181) q[2];
rz(1.5859531) q[3];
sx q[3];
rz(-2.0944244) q[3];
sx q[3];
rz(-1.9716523) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8856119) q[0];
sx q[0];
rz(-2.8035127) q[0];
sx q[0];
rz(1.0611435) q[0];
rz(1.9288918) q[1];
sx q[1];
rz(-2.3586528) q[1];
sx q[1];
rz(-1.3139668) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6550094) q[0];
sx q[0];
rz(-0.81196852) q[0];
sx q[0];
rz(2.6467108) q[0];
rz(-2.2167614) q[2];
sx q[2];
rz(-1.2154757) q[2];
sx q[2];
rz(-2.3424847) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3281341) q[1];
sx q[1];
rz(-2.5812006) q[1];
sx q[1];
rz(1.2499534) q[1];
rz(0.84557326) q[3];
sx q[3];
rz(-1.9872287) q[3];
sx q[3];
rz(1.6628979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5709915) q[2];
sx q[2];
rz(-0.97390318) q[2];
sx q[2];
rz(-0.904733) q[2];
rz(1.5915126) q[3];
sx q[3];
rz(-1.2856057) q[3];
sx q[3];
rz(1.2929644) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631113) q[0];
sx q[0];
rz(-3.0516629) q[0];
sx q[0];
rz(2.2233546) q[0];
rz(-0.69084424) q[1];
sx q[1];
rz(-1.3965239) q[1];
sx q[1];
rz(0.7965368) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.349571) q[0];
sx q[0];
rz(-2.9137028) q[0];
sx q[0];
rz(-1.9877276) q[0];
rz(0.94743016) q[2];
sx q[2];
rz(-0.073270144) q[2];
sx q[2];
rz(2.2712592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0738363) q[1];
sx q[1];
rz(-0.87693611) q[1];
sx q[1];
rz(0.14090726) q[1];
rz(-2.2240766) q[3];
sx q[3];
rz(-1.8492572) q[3];
sx q[3];
rz(2.3115932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.010926509) q[2];
sx q[2];
rz(-0.9504168) q[2];
sx q[2];
rz(1.2476904) q[2];
rz(-0.55825663) q[3];
sx q[3];
rz(-1.6853251) q[3];
sx q[3];
rz(-3.1202417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.438544) q[0];
sx q[0];
rz(-3.092364) q[0];
sx q[0];
rz(-2.4141648) q[0];
rz(-2.9797331) q[1];
sx q[1];
rz(-1.1715803) q[1];
sx q[1];
rz(-1.1579827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25140554) q[0];
sx q[0];
rz(-1.5047538) q[0];
sx q[0];
rz(1.4736045) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.00096614758) q[2];
sx q[2];
rz(-1.4314993) q[2];
sx q[2];
rz(-1.7666221) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.53044415) q[1];
sx q[1];
rz(-1.6940677) q[1];
sx q[1];
rz(1.7632496) q[1];
rz(-pi) q[2];
rz(-0.14258464) q[3];
sx q[3];
rz(-0.53053427) q[3];
sx q[3];
rz(-0.50850463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8689279) q[2];
sx q[2];
rz(-1.9395892) q[2];
sx q[2];
rz(1.7499917) q[2];
rz(0.23022716) q[3];
sx q[3];
rz(-1.1743436) q[3];
sx q[3];
rz(-0.19190425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72628438) q[0];
sx q[0];
rz(-0.10157651) q[0];
sx q[0];
rz(-2.0368077) q[0];
rz(-0.11100189) q[1];
sx q[1];
rz(-1.1155201) q[1];
sx q[1];
rz(1.312779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3462145) q[0];
sx q[0];
rz(-1.7722478) q[0];
sx q[0];
rz(0.18152118) q[0];
rz(-0.13938015) q[2];
sx q[2];
rz(-0.5682033) q[2];
sx q[2];
rz(-2.2903493) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5178419) q[1];
sx q[1];
rz(-1.1718796) q[1];
sx q[1];
rz(2.9241882) q[1];
rz(-pi) q[2];
rz(-0.41683414) q[3];
sx q[3];
rz(-2.2456944) q[3];
sx q[3];
rz(1.0639497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63867265) q[2];
sx q[2];
rz(-1.0588249) q[2];
sx q[2];
rz(2.1144833) q[2];
rz(-2.1549759) q[3];
sx q[3];
rz(-0.28237453) q[3];
sx q[3];
rz(-1.4237039) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8091549) q[0];
sx q[0];
rz(-1.2323392) q[0];
sx q[0];
rz(-0.80068457) q[0];
rz(-2.370131) q[1];
sx q[1];
rz(-1.0779251) q[1];
sx q[1];
rz(-1.7652184) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3478394) q[0];
sx q[0];
rz(-2.3476566) q[0];
sx q[0];
rz(-0.57200045) q[0];
rz(-2.0368783) q[2];
sx q[2];
rz(-1.5395293) q[2];
sx q[2];
rz(-0.75541562) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.034198) q[1];
sx q[1];
rz(-1.5840816) q[1];
sx q[1];
rz(-0.70550578) q[1];
rz(2.3153051) q[3];
sx q[3];
rz(-2.1231696) q[3];
sx q[3];
rz(-0.9780405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2996404) q[2];
sx q[2];
rz(-1.2066634) q[2];
sx q[2];
rz(-0.017814962) q[2];
rz(-0.31339112) q[3];
sx q[3];
rz(-1.0217977) q[3];
sx q[3];
rz(2.4221086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3571091) q[0];
sx q[0];
rz(-1.5743558) q[0];
sx q[0];
rz(2.9743279) q[0];
rz(-0.74370614) q[1];
sx q[1];
rz(-0.88631648) q[1];
sx q[1];
rz(3.0053265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0038827) q[0];
sx q[0];
rz(-1.2943177) q[0];
sx q[0];
rz(2.2599758) q[0];
x q[1];
rz(2.9685814) q[2];
sx q[2];
rz(-1.271651) q[2];
sx q[2];
rz(-2.0188221) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4748342) q[1];
sx q[1];
rz(-1.5038361) q[1];
sx q[1];
rz(-2.7826742) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6266009) q[3];
sx q[3];
rz(-1.9161092) q[3];
sx q[3];
rz(1.8453516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83884376) q[2];
sx q[2];
rz(-0.15041298) q[2];
sx q[2];
rz(2.0602843) q[2];
rz(-2.9980764) q[3];
sx q[3];
rz(-1.84294) q[3];
sx q[3];
rz(1.4474086) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3926587) q[0];
sx q[0];
rz(-1.3329788) q[0];
sx q[0];
rz(-0.64312154) q[0];
rz(-2.2195623) q[1];
sx q[1];
rz(-2.8755201) q[1];
sx q[1];
rz(1.5111074) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9859146) q[0];
sx q[0];
rz(-1.2307457) q[0];
sx q[0];
rz(-2.9513533) q[0];
rz(1.8505356) q[2];
sx q[2];
rz(-1.4609173) q[2];
sx q[2];
rz(2.2834275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0865759) q[1];
sx q[1];
rz(-0.19682238) q[1];
sx q[1];
rz(1.5301401) q[1];
x q[2];
rz(1.891252) q[3];
sx q[3];
rz(-2.0600187) q[3];
sx q[3];
rz(1.7699522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1954605) q[2];
sx q[2];
rz(-1.5798502) q[2];
sx q[2];
rz(-1.8411609) q[2];
rz(-0.91222936) q[3];
sx q[3];
rz(-1.8650841) q[3];
sx q[3];
rz(-2.8279878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6533971) q[0];
sx q[0];
rz(-2.6618239) q[0];
sx q[0];
rz(-2.524014) q[0];
rz(-3.0600582) q[1];
sx q[1];
rz(-0.50059861) q[1];
sx q[1];
rz(-2.4542782) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9195557) q[0];
sx q[0];
rz(-1.3341781) q[0];
sx q[0];
rz(0.048892269) q[0];
rz(1.4590644) q[2];
sx q[2];
rz(-1.8364753) q[2];
sx q[2];
rz(0.39171644) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3442823) q[1];
sx q[1];
rz(-0.51206368) q[1];
sx q[1];
rz(-0.6060098) q[1];
rz(-3.0932759) q[3];
sx q[3];
rz(-2.1135277) q[3];
sx q[3];
rz(0.84191546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0240137) q[2];
sx q[2];
rz(-1.964183) q[2];
sx q[2];
rz(3.0702316) q[2];
rz(1.2355545) q[3];
sx q[3];
rz(-0.33696285) q[3];
sx q[3];
rz(-2.5753042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5736893) q[0];
sx q[0];
rz(-2.8901143) q[0];
sx q[0];
rz(-0.13791826) q[0];
rz(0.57016405) q[1];
sx q[1];
rz(-2.0591996) q[1];
sx q[1];
rz(-0.015448419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.100525) q[0];
sx q[0];
rz(-2.3599632) q[0];
sx q[0];
rz(1.1055205) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3809468) q[2];
sx q[2];
rz(-1.9234675) q[2];
sx q[2];
rz(-1.901153) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8739104) q[1];
sx q[1];
rz(-1.7355669) q[1];
sx q[1];
rz(2.6450637) q[1];
rz(-pi) q[2];
rz(0.055850765) q[3];
sx q[3];
rz(-1.3875563) q[3];
sx q[3];
rz(2.2574772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7323759) q[2];
sx q[2];
rz(-2.3541383) q[2];
sx q[2];
rz(-2.709205) q[2];
rz(0.89808291) q[3];
sx q[3];
rz(-0.87528527) q[3];
sx q[3];
rz(1.557365) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6569923) q[0];
sx q[0];
rz(-2.0130172) q[0];
sx q[0];
rz(0.76242557) q[0];
rz(-2.5905329) q[1];
sx q[1];
rz(-1.8012128) q[1];
sx q[1];
rz(2.6691379) q[1];
rz(-3.0226784) q[2];
sx q[2];
rz(-0.36570315) q[2];
sx q[2];
rz(-0.23592532) q[2];
rz(-0.9394218) q[3];
sx q[3];
rz(-0.71790725) q[3];
sx q[3];
rz(0.090286615) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
