OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8853814) q[0];
sx q[0];
rz(-2.4282832) q[0];
sx q[0];
rz(-1.927884) q[0];
rz(-1.3656536) q[1];
sx q[1];
rz(-2.8627099) q[1];
sx q[1];
rz(-0.20067781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047223481) q[0];
sx q[0];
rz(-2.6728164) q[0];
sx q[0];
rz(-1.7550521) q[0];
rz(-1.7332478) q[2];
sx q[2];
rz(-1.7853569) q[2];
sx q[2];
rz(-1.868737) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9057718) q[1];
sx q[1];
rz(-1.3911472) q[1];
sx q[1];
rz(0.19975042) q[1];
rz(-pi) q[2];
rz(2.8332769) q[3];
sx q[3];
rz(-1.7349958) q[3];
sx q[3];
rz(-0.37946821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.063244907) q[2];
sx q[2];
rz(-1.1996317) q[2];
sx q[2];
rz(-0.71301785) q[2];
rz(-1.2837422) q[3];
sx q[3];
rz(-0.72834891) q[3];
sx q[3];
rz(-0.89948765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9707608) q[0];
sx q[0];
rz(-1.4749227) q[0];
sx q[0];
rz(-1.4738039) q[0];
rz(2.5693192) q[1];
sx q[1];
rz(-2.0794561) q[1];
sx q[1];
rz(-1.7639814) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77967465) q[0];
sx q[0];
rz(-0.55680823) q[0];
sx q[0];
rz(-0.98018719) q[0];
x q[1];
rz(0.53133659) q[2];
sx q[2];
rz(-1.9730933) q[2];
sx q[2];
rz(0.50237331) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2050537) q[1];
sx q[1];
rz(-1.8636101) q[1];
sx q[1];
rz(1.9163301) q[1];
x q[2];
rz(-0.20444025) q[3];
sx q[3];
rz(-1.2693161) q[3];
sx q[3];
rz(0.35273509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1231692) q[2];
sx q[2];
rz(-1.5095242) q[2];
sx q[2];
rz(1.8040166) q[2];
rz(0.34559524) q[3];
sx q[3];
rz(-2.5497422) q[3];
sx q[3];
rz(0.93446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086394101) q[0];
sx q[0];
rz(-2.9559657) q[0];
sx q[0];
rz(0.58919543) q[0];
rz(-0.20801726) q[1];
sx q[1];
rz(-2.5841525) q[1];
sx q[1];
rz(-0.20588188) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1070646) q[0];
sx q[0];
rz(-2.42332) q[0];
sx q[0];
rz(-2.5884401) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0127064) q[2];
sx q[2];
rz(-2.3461902) q[2];
sx q[2];
rz(-2.263042) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.65218788) q[1];
sx q[1];
rz(-0.88880782) q[1];
sx q[1];
rz(-1.9661862) q[1];
rz(1.7725189) q[3];
sx q[3];
rz(-0.414114) q[3];
sx q[3];
rz(2.5270346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0516633) q[2];
sx q[2];
rz(-1.0445003) q[2];
sx q[2];
rz(-0.54452407) q[2];
rz(0.45474592) q[3];
sx q[3];
rz(-0.62362042) q[3];
sx q[3];
rz(-0.0039984306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74314463) q[0];
sx q[0];
rz(-0.54773098) q[0];
sx q[0];
rz(2.4416583) q[0];
rz(-1.8327389) q[1];
sx q[1];
rz(-2.0715641) q[1];
sx q[1];
rz(1.4189789) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0467906) q[0];
sx q[0];
rz(-1.3700458) q[0];
sx q[0];
rz(2.0457474) q[0];
x q[1];
rz(1.5350902) q[2];
sx q[2];
rz(-2.6427445) q[2];
sx q[2];
rz(-0.6998261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0154889) q[1];
sx q[1];
rz(-0.88147336) q[1];
sx q[1];
rz(1.8828805) q[1];
rz(2.9955825) q[3];
sx q[3];
rz(-1.8838804) q[3];
sx q[3];
rz(-2.2816471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2122638) q[2];
sx q[2];
rz(-2.6876891) q[2];
sx q[2];
rz(1.3455428) q[2];
rz(0.71873194) q[3];
sx q[3];
rz(-2.4001382) q[3];
sx q[3];
rz(-0.68812266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6253925) q[0];
sx q[0];
rz(-1.9560408) q[0];
sx q[0];
rz(2.5984883) q[0];
rz(-0.25554666) q[1];
sx q[1];
rz(-2.6592022) q[1];
sx q[1];
rz(1.3822752) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2204748) q[0];
sx q[0];
rz(-0.68604031) q[0];
sx q[0];
rz(-0.6612079) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9901401) q[2];
sx q[2];
rz(-2.3690358) q[2];
sx q[2];
rz(-2.549215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2321736) q[1];
sx q[1];
rz(-0.62271732) q[1];
sx q[1];
rz(-1.1145341) q[1];
rz(-1.3997283) q[3];
sx q[3];
rz(-1.1753963) q[3];
sx q[3];
rz(-2.4912437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1725267) q[2];
sx q[2];
rz(-1.4719937) q[2];
sx q[2];
rz(-2.772061) q[2];
rz(-1.3299804) q[3];
sx q[3];
rz(-2.3510635) q[3];
sx q[3];
rz(-1.3548405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0238277) q[0];
sx q[0];
rz(-2.8067639) q[0];
sx q[0];
rz(-3.0552926) q[0];
rz(1.6465126) q[1];
sx q[1];
rz(-1.1777271) q[1];
sx q[1];
rz(0.42974791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1275009) q[0];
sx q[0];
rz(-1.3821116) q[0];
sx q[0];
rz(-3.0701915) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8977857) q[2];
sx q[2];
rz(-2.2026416) q[2];
sx q[2];
rz(-1.6666842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60372231) q[1];
sx q[1];
rz(-0.66927823) q[1];
sx q[1];
rz(-0.75103384) q[1];
x q[2];
rz(0.34853156) q[3];
sx q[3];
rz(-2.4405839) q[3];
sx q[3];
rz(1.9402869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0325844) q[2];
sx q[2];
rz(-2.0592368) q[2];
sx q[2];
rz(-1.9411795) q[2];
rz(2.9754908) q[3];
sx q[3];
rz(-0.76986543) q[3];
sx q[3];
rz(2.2782245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2849041) q[0];
sx q[0];
rz(-2.3785474) q[0];
sx q[0];
rz(-0.82409182) q[0];
rz(1.2245945) q[1];
sx q[1];
rz(-1.2242182) q[1];
sx q[1];
rz(1.0138938) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4201502) q[0];
sx q[0];
rz(-2.3334896) q[0];
sx q[0];
rz(1.1366913) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8877007) q[2];
sx q[2];
rz(-1.1528258) q[2];
sx q[2];
rz(2.5078513) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6142004) q[1];
sx q[1];
rz(-2.1062615) q[1];
sx q[1];
rz(-2.0245069) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1070293) q[3];
sx q[3];
rz(-1.9120029) q[3];
sx q[3];
rz(-0.1078913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8635204) q[2];
sx q[2];
rz(-0.72276989) q[2];
sx q[2];
rz(1.8966759) q[2];
rz(-0.23166367) q[3];
sx q[3];
rz(-2.1004227) q[3];
sx q[3];
rz(-2.9932573) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0949377) q[0];
sx q[0];
rz(-1.2060839) q[0];
sx q[0];
rz(-0.99223247) q[0];
rz(0.16920432) q[1];
sx q[1];
rz(-2.2997586) q[1];
sx q[1];
rz(1.7281035) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9887392) q[0];
sx q[0];
rz(-2.1655575) q[0];
sx q[0];
rz(1.4805984) q[0];
rz(-pi) q[1];
rz(-2.1229073) q[2];
sx q[2];
rz(-1.3231965) q[2];
sx q[2];
rz(-1.1820861) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38138776) q[1];
sx q[1];
rz(-1.5389812) q[1];
sx q[1];
rz(2.8434668) q[1];
rz(-2.9417073) q[3];
sx q[3];
rz(-2.0428786) q[3];
sx q[3];
rz(-0.98552824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8168872) q[2];
sx q[2];
rz(-1.6605261) q[2];
sx q[2];
rz(1.7186349) q[2];
rz(0.11296806) q[3];
sx q[3];
rz(-0.26765099) q[3];
sx q[3];
rz(0.62937361) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9519575) q[0];
sx q[0];
rz(-1.676214) q[0];
sx q[0];
rz(-2.5919609) q[0];
rz(2.5426087) q[1];
sx q[1];
rz(-2.2689029) q[1];
sx q[1];
rz(1.6302861) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0534759) q[0];
sx q[0];
rz(-1.5742745) q[0];
sx q[0];
rz(-0.17668488) q[0];
rz(-2.6341166) q[2];
sx q[2];
rz(-1.8364454) q[2];
sx q[2];
rz(-1.332559) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0508742) q[1];
sx q[1];
rz(-1.1696739) q[1];
sx q[1];
rz(-2.119333) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6502147) q[3];
sx q[3];
rz(-1.8432321) q[3];
sx q[3];
rz(-2.2898554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1134243) q[2];
sx q[2];
rz(-1.7999962) q[2];
sx q[2];
rz(2.8070731) q[2];
rz(2.4962375) q[3];
sx q[3];
rz(-1.0048451) q[3];
sx q[3];
rz(0.2373124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0765814) q[0];
sx q[0];
rz(-0.33877057) q[0];
sx q[0];
rz(0.6231128) q[0];
rz(-0.30162853) q[1];
sx q[1];
rz(-0.61768618) q[1];
sx q[1];
rz(2.1004486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17693263) q[0];
sx q[0];
rz(-0.76632351) q[0];
sx q[0];
rz(0.64789741) q[0];
rz(-pi) q[1];
rz(-0.1847965) q[2];
sx q[2];
rz(-2.6738538) q[2];
sx q[2];
rz(0.50450215) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5981982) q[1];
sx q[1];
rz(-2.3375727) q[1];
sx q[1];
rz(-1.6135741) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9277293) q[3];
sx q[3];
rz(-1.2215316) q[3];
sx q[3];
rz(1.6598778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4447896) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(-0.94318548) q[2];
rz(1.7275564) q[3];
sx q[3];
rz(-0.90641886) q[3];
sx q[3];
rz(2.3448155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74054756) q[0];
sx q[0];
rz(-0.39401207) q[0];
sx q[0];
rz(3.0643585) q[0];
rz(-0.41863353) q[1];
sx q[1];
rz(-1.5593465) q[1];
sx q[1];
rz(-2.9895463) q[1];
rz(1.536191) q[2];
sx q[2];
rz(-0.99609756) q[2];
sx q[2];
rz(0.20142557) q[2];
rz(-0.54316212) q[3];
sx q[3];
rz(-0.49351963) q[3];
sx q[3];
rz(-1.7388572) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
