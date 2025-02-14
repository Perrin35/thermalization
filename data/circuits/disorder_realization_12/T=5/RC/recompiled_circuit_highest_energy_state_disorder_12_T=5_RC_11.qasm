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
rz(-0.15859088) q[0];
sx q[0];
rz(-0.51704419) q[0];
sx q[0];
rz(1.3102732) q[0];
rz(2.7994371) q[1];
sx q[1];
rz(-1.9603536) q[1];
sx q[1];
rz(-0.96460834) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5072381) q[0];
sx q[0];
rz(-0.94598457) q[0];
sx q[0];
rz(-1.4039197) q[0];
rz(-pi) q[1];
rz(1.3753424) q[2];
sx q[2];
rz(-1.9103622) q[2];
sx q[2];
rz(-1.5887141) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1312932) q[1];
sx q[1];
rz(-2.2939766) q[1];
sx q[1];
rz(2.7406969) q[1];
x q[2];
rz(0.11537376) q[3];
sx q[3];
rz(-0.90297195) q[3];
sx q[3];
rz(-1.7911719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0059119314) q[2];
sx q[2];
rz(-2.9475806) q[2];
sx q[2];
rz(1.4646336) q[2];
rz(-1.3095193) q[3];
sx q[3];
rz(-0.94683164) q[3];
sx q[3];
rz(-0.32002282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86318535) q[0];
sx q[0];
rz(-1.139737) q[0];
sx q[0];
rz(1.3673258) q[0];
rz(-1.6525846) q[1];
sx q[1];
rz(-2.3431578) q[1];
sx q[1];
rz(2.8489825) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6824043) q[0];
sx q[0];
rz(-1.4727797) q[0];
sx q[0];
rz(-2.8062264) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9797249) q[2];
sx q[2];
rz(-1.6965995) q[2];
sx q[2];
rz(2.6243072) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78116592) q[1];
sx q[1];
rz(-1.806621) q[1];
sx q[1];
rz(1.3597914) q[1];
rz(2.9037016) q[3];
sx q[3];
rz(-2.8679016) q[3];
sx q[3];
rz(-2.9546933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9802398) q[2];
sx q[2];
rz(-0.62675256) q[2];
sx q[2];
rz(-2.9577067) q[2];
rz(0.45267496) q[3];
sx q[3];
rz(-1.6190745) q[3];
sx q[3];
rz(3.0116853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1995131) q[0];
sx q[0];
rz(-2.4266854) q[0];
sx q[0];
rz(2.3364501) q[0];
rz(2.0623656) q[1];
sx q[1];
rz(-2.3715623) q[1];
sx q[1];
rz(-1.8260746) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0502687) q[0];
sx q[0];
rz(-1.5116232) q[0];
sx q[0];
rz(2.4135804) q[0];
rz(-2.3147624) q[2];
sx q[2];
rz(-1.9212124) q[2];
sx q[2];
rz(1.0817127) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6650369) q[1];
sx q[1];
rz(-1.2893234) q[1];
sx q[1];
rz(-3.046077) q[1];
x q[2];
rz(1.4844378) q[3];
sx q[3];
rz(-1.8724073) q[3];
sx q[3];
rz(-2.3197966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2074073) q[2];
sx q[2];
rz(-1.6560053) q[2];
sx q[2];
rz(-0.54971131) q[2];
rz(-0.011119757) q[3];
sx q[3];
rz(-2.34237) q[3];
sx q[3];
rz(-1.3219249) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84871197) q[0];
sx q[0];
rz(-0.95893812) q[0];
sx q[0];
rz(-2.8380561) q[0];
rz(0.15829463) q[1];
sx q[1];
rz(-2.0469432) q[1];
sx q[1];
rz(1.0096445) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15561129) q[0];
sx q[0];
rz(-1.2878875) q[0];
sx q[0];
rz(2.020218) q[0];
x q[1];
rz(2.7042806) q[2];
sx q[2];
rz(-0.6491937) q[2];
sx q[2];
rz(-1.5199666) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18100723) q[1];
sx q[1];
rz(-0.50324856) q[1];
sx q[1];
rz(-2.8099847) q[1];
rz(-0.03831717) q[3];
sx q[3];
rz(-2.6434757) q[3];
sx q[3];
rz(1.8565053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2046854) q[2];
sx q[2];
rz(-2.3202809) q[2];
sx q[2];
rz(-0.26910195) q[2];
rz(1.4540295) q[3];
sx q[3];
rz(-1.5498091) q[3];
sx q[3];
rz(1.5788186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8708385) q[0];
sx q[0];
rz(-1.0601059) q[0];
sx q[0];
rz(0.36218542) q[0];
rz(-2.4902792) q[1];
sx q[1];
rz(-1.6505227) q[1];
sx q[1];
rz(1.5247033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3782048) q[0];
sx q[0];
rz(-1.7001061) q[0];
sx q[0];
rz(2.5300701) q[0];
rz(-pi) q[1];
rz(-0.082090898) q[2];
sx q[2];
rz(-0.94952784) q[2];
sx q[2];
rz(0.7393078) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0901047) q[1];
sx q[1];
rz(-1.0516775) q[1];
sx q[1];
rz(0.41579065) q[1];
x q[2];
rz(2.7157341) q[3];
sx q[3];
rz(-0.79090624) q[3];
sx q[3];
rz(-0.87825852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1093971) q[2];
sx q[2];
rz(-0.95462644) q[2];
sx q[2];
rz(-1.7677914) q[2];
rz(2.7023756) q[3];
sx q[3];
rz(-2.4924811) q[3];
sx q[3];
rz(-2.5325328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.1010901) q[0];
sx q[0];
rz(-2.4093565) q[0];
sx q[0];
rz(2.8771583) q[0];
rz(0.6791555) q[1];
sx q[1];
rz(-0.87398386) q[1];
sx q[1];
rz(0.98495475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20709461) q[0];
sx q[0];
rz(-2.6885928) q[0];
sx q[0];
rz(-2.3191602) q[0];
rz(2.856065) q[2];
sx q[2];
rz(-2.0047024) q[2];
sx q[2];
rz(1.3615695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90675844) q[1];
sx q[1];
rz(-1.088716) q[1];
sx q[1];
rz(1.158403) q[1];
rz(-pi) q[2];
rz(1.1285176) q[3];
sx q[3];
rz(-1.3849713) q[3];
sx q[3];
rz(3.1215661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1017477) q[2];
sx q[2];
rz(-0.16182772) q[2];
sx q[2];
rz(0.43231535) q[2];
rz(1.013422) q[3];
sx q[3];
rz(-1.5347967) q[3];
sx q[3];
rz(0.54949808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48563114) q[0];
sx q[0];
rz(-0.39325842) q[0];
sx q[0];
rz(2.0320758) q[0];
rz(-2.3484777) q[1];
sx q[1];
rz(-1.3536645) q[1];
sx q[1];
rz(1.5464334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30013672) q[0];
sx q[0];
rz(-1.0086035) q[0];
sx q[0];
rz(-2.4676222) q[0];
rz(1.7810473) q[2];
sx q[2];
rz(-0.86019197) q[2];
sx q[2];
rz(0.17786121) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2364051) q[1];
sx q[1];
rz(-1.3705472) q[1];
sx q[1];
rz(-0.88242857) q[1];
x q[2];
rz(-0.1216573) q[3];
sx q[3];
rz(-0.72008789) q[3];
sx q[3];
rz(0.64313221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3194797) q[2];
sx q[2];
rz(-2.959368) q[2];
sx q[2];
rz(1.7721843) q[2];
rz(-2.2111514) q[3];
sx q[3];
rz(-1.3481827) q[3];
sx q[3];
rz(1.6377009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9419788) q[0];
sx q[0];
rz(-1.9604585) q[0];
sx q[0];
rz(-0.30366316) q[0];
rz(2.1619201) q[1];
sx q[1];
rz(-0.62398463) q[1];
sx q[1];
rz(0.40506515) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41202085) q[0];
sx q[0];
rz(-1.7639909) q[0];
sx q[0];
rz(-2.7576819) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.151565) q[2];
sx q[2];
rz(-1.4128601) q[2];
sx q[2];
rz(-1.7309703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0593589) q[1];
sx q[1];
rz(-0.59525604) q[1];
sx q[1];
rz(2.8221376) q[1];
x q[2];
rz(-1.1438683) q[3];
sx q[3];
rz(-1.3865304) q[3];
sx q[3];
rz(-1.9686521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27641174) q[2];
sx q[2];
rz(-0.96651912) q[2];
sx q[2];
rz(-2.442404) q[2];
rz(-0.80896038) q[3];
sx q[3];
rz(-1.304108) q[3];
sx q[3];
rz(0.32304025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9998099) q[0];
sx q[0];
rz(-1.9843822) q[0];
sx q[0];
rz(-3.004177) q[0];
rz(-3.0925062) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(2.633599) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2823454) q[0];
sx q[0];
rz(-2.0657592) q[0];
sx q[0];
rz(1.540394) q[0];
x q[1];
rz(-1.5586583) q[2];
sx q[2];
rz(-1.3985139) q[2];
sx q[2];
rz(-2.7207295) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55785364) q[1];
sx q[1];
rz(-2.5370829) q[1];
sx q[1];
rz(2.1624751) q[1];
rz(3.052086) q[3];
sx q[3];
rz(-1.4495159) q[3];
sx q[3];
rz(-2.3119777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55234838) q[2];
sx q[2];
rz(-1.9231223) q[2];
sx q[2];
rz(-2.153896) q[2];
rz(0.47932953) q[3];
sx q[3];
rz(-1.597581) q[3];
sx q[3];
rz(-2.2492669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
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
rz(2.4722897) q[0];
sx q[0];
rz(-0.23858128) q[0];
sx q[0];
rz(0.14078374) q[0];
rz(-0.36551481) q[1];
sx q[1];
rz(-2.3194158) q[1];
sx q[1];
rz(-0.64868322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1295531) q[0];
sx q[0];
rz(-2.3863992) q[0];
sx q[0];
rz(-2.6213403) q[0];
rz(-pi) q[1];
rz(-3.0200483) q[2];
sx q[2];
rz(-0.76833188) q[2];
sx q[2];
rz(-2.2939081) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4280871) q[1];
sx q[1];
rz(-2.1660342) q[1];
sx q[1];
rz(-1.6387322) q[1];
rz(-2.8148531) q[3];
sx q[3];
rz(-1.3413197) q[3];
sx q[3];
rz(-0.47577259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6792128) q[2];
sx q[2];
rz(-1.0054192) q[2];
sx q[2];
rz(-3.0688378) q[2];
rz(0.66271979) q[3];
sx q[3];
rz(-1.8279671) q[3];
sx q[3];
rz(-0.16547468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0543906) q[0];
sx q[0];
rz(-1.6086171) q[0];
sx q[0];
rz(1.4679012) q[0];
rz(1.5351334) q[1];
sx q[1];
rz(-0.88773334) q[1];
sx q[1];
rz(1.6233374) q[1];
rz(0.4735422) q[2];
sx q[2];
rz(-2.8859856) q[2];
sx q[2];
rz(-2.4398266) q[2];
rz(0.88964626) q[3];
sx q[3];
rz(-1.2500661) q[3];
sx q[3];
rz(1.4483857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
