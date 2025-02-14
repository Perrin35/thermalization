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
rz(1.5155091) q[0];
sx q[0];
rz(3.4184366) q[0];
sx q[0];
rz(9.3079216) q[0];
rz(2.6729743) q[1];
sx q[1];
rz(-0.35419551) q[1];
sx q[1];
rz(0.49122223) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15165937) q[0];
sx q[0];
rz(-1.0826246) q[0];
sx q[0];
rz(-3.0753646) q[0];
rz(2.6543219) q[2];
sx q[2];
rz(-1.7585764) q[2];
sx q[2];
rz(-2.1405381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3892606) q[1];
sx q[1];
rz(-1.4717143) q[1];
sx q[1];
rz(-2.5072751) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9141044) q[3];
sx q[3];
rz(-1.5945537) q[3];
sx q[3];
rz(-2.6809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7626071) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(0.22988764) q[2];
rz(3.0843132) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(0.91276401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7951935) q[0];
sx q[0];
rz(-2.754358) q[0];
sx q[0];
rz(-0.21695319) q[0];
rz(3.1087061) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(-2.4864973) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0467619) q[0];
sx q[0];
rz(-0.74709409) q[0];
sx q[0];
rz(-1.2741035) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7704324) q[2];
sx q[2];
rz(-2.0642274) q[2];
sx q[2];
rz(2.9889929) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2456783) q[1];
sx q[1];
rz(-2.0729471) q[1];
sx q[1];
rz(0.72523004) q[1];
x q[2];
rz(2.3052196) q[3];
sx q[3];
rz(-1.4640088) q[3];
sx q[3];
rz(1.2272782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32404262) q[2];
sx q[2];
rz(-2.5174759) q[2];
sx q[2];
rz(-1.2980655) q[2];
rz(-1.2201747) q[3];
sx q[3];
rz(-1.7036567) q[3];
sx q[3];
rz(-0.65742457) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87797457) q[0];
sx q[0];
rz(-0.89732301) q[0];
sx q[0];
rz(-2.5275912) q[0];
rz(1.7738495) q[1];
sx q[1];
rz(-2.4847993) q[1];
sx q[1];
rz(0.49555379) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14392631) q[0];
sx q[0];
rz(-2.196402) q[0];
sx q[0];
rz(2.1493874) q[0];
rz(-pi) q[1];
rz(2.6265385) q[2];
sx q[2];
rz(-0.94991377) q[2];
sx q[2];
rz(-2.1504998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6537955) q[1];
sx q[1];
rz(-1.3843517) q[1];
sx q[1];
rz(1.3273147) q[1];
rz(0.46681459) q[3];
sx q[3];
rz(-1.8895849) q[3];
sx q[3];
rz(2.7960645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76501781) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(-0.59256727) q[2];
rz(2.9060034) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(-2.6872915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6465103) q[0];
sx q[0];
rz(-0.80146924) q[0];
sx q[0];
rz(0.0058280514) q[0];
rz(0.5799154) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(2.6204806) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2297008) q[0];
sx q[0];
rz(-1.5754412) q[0];
sx q[0];
rz(1.5779431) q[0];
rz(-2.1717511) q[2];
sx q[2];
rz(-1.5464939) q[2];
sx q[2];
rz(1.5497335) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5856984) q[1];
sx q[1];
rz(-1.1194849) q[1];
sx q[1];
rz(2.0183619) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24346015) q[3];
sx q[3];
rz(-0.32103466) q[3];
sx q[3];
rz(2.3784172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8899322) q[2];
sx q[2];
rz(-2.6106788) q[2];
sx q[2];
rz(-0.31822515) q[2];
rz(-2.4288154) q[3];
sx q[3];
rz(-0.30064279) q[3];
sx q[3];
rz(1.4819063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37079674) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(0.036291432) q[0];
rz(2.6151784) q[1];
sx q[1];
rz(-1.4585835) q[1];
sx q[1];
rz(0.2821736) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.605114) q[0];
sx q[0];
rz(-1.295034) q[0];
sx q[0];
rz(0.92425297) q[0];
x q[1];
rz(2.2660115) q[2];
sx q[2];
rz(-1.1291847) q[2];
sx q[2];
rz(1.832806) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1212738) q[1];
sx q[1];
rz(-1.6409546) q[1];
sx q[1];
rz(2.4362045) q[1];
x q[2];
rz(-2.3822278) q[3];
sx q[3];
rz(-2.1753484) q[3];
sx q[3];
rz(0.60455632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37461773) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(0.10147632) q[2];
rz(0.96595079) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(-0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517956) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(1.1740603) q[0];
rz(-0.38408285) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(-0.039965872) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46972457) q[0];
sx q[0];
rz(-1.3626119) q[0];
sx q[0];
rz(-1.2428268) q[0];
x q[1];
rz(2.2117679) q[2];
sx q[2];
rz(-1.3017968) q[2];
sx q[2];
rz(-2.9421634) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.080885305) q[1];
sx q[1];
rz(-0.97700143) q[1];
sx q[1];
rz(1.9869366) q[1];
x q[2];
rz(-1.3371991) q[3];
sx q[3];
rz(-0.75832483) q[3];
sx q[3];
rz(2.9616464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93969718) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(0.8197909) q[2];
rz(0.78153265) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(1.3801105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.49067295) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(-1.3167205) q[0];
rz(-1.7735749) q[1];
sx q[1];
rz(-1.2317692) q[1];
sx q[1];
rz(-2.7046611) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9428944) q[0];
sx q[0];
rz(-1.534053) q[0];
sx q[0];
rz(-1.3803287) q[0];
rz(-pi) q[1];
rz(-1.6810136) q[2];
sx q[2];
rz(-1.6880562) q[2];
sx q[2];
rz(0.77484581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93561427) q[1];
sx q[1];
rz(-0.84594489) q[1];
sx q[1];
rz(-2.0589252) q[1];
rz(-0.42610069) q[3];
sx q[3];
rz(-1.0747194) q[3];
sx q[3];
rz(-2.8471198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4009319) q[2];
sx q[2];
rz(-2.0621982) q[2];
sx q[2];
rz(0.13460049) q[2];
rz(1.918119) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(1.9756165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3400035) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(-0.18184161) q[0];
rz(0.34135154) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(-2.9827319) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7094966) q[0];
sx q[0];
rz(-2.0176945) q[0];
sx q[0];
rz(1.4947557) q[0];
rz(-2.5652021) q[2];
sx q[2];
rz(-0.49242556) q[2];
sx q[2];
rz(-0.8053636) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78217172) q[1];
sx q[1];
rz(-1.8295994) q[1];
sx q[1];
rz(1.7951629) q[1];
rz(-pi) q[2];
rz(-2.587593) q[3];
sx q[3];
rz(-1.7007728) q[3];
sx q[3];
rz(0.25773559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2813256) q[2];
sx q[2];
rz(-0.93026668) q[2];
sx q[2];
rz(-2.2567828) q[2];
rz(-1.7224711) q[3];
sx q[3];
rz(-2.116674) q[3];
sx q[3];
rz(2.6917246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.849843) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(-2.1480609) q[0];
rz(-1.7811071) q[1];
sx q[1];
rz(-2.2527756) q[1];
sx q[1];
rz(0.55307585) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4489514) q[0];
sx q[0];
rz(-1.7990489) q[0];
sx q[0];
rz(-0.3897764) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73846784) q[2];
sx q[2];
rz(-1.6203801) q[2];
sx q[2];
rz(2.2253583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.304405) q[1];
sx q[1];
rz(-2.4585729) q[1];
sx q[1];
rz(-2.7262615) q[1];
x q[2];
rz(-0.27121765) q[3];
sx q[3];
rz(-2.1925049) q[3];
sx q[3];
rz(0.78628899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9370344) q[2];
sx q[2];
rz(-0.60595804) q[2];
sx q[2];
rz(0.23703144) q[2];
rz(-1.5406746) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(-1.1698394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0700664) q[0];
sx q[0];
rz(-1.5109753) q[0];
sx q[0];
rz(3.0226829) q[0];
rz(-0.43926829) q[1];
sx q[1];
rz(-1.2989137) q[1];
sx q[1];
rz(3.0781504) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251415) q[0];
sx q[0];
rz(-2.2277895) q[0];
sx q[0];
rz(-2.444066) q[0];
rz(-pi) q[1];
rz(0.25524885) q[2];
sx q[2];
rz(-2.5016157) q[2];
sx q[2];
rz(0.6682804) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0234785) q[1];
sx q[1];
rz(-2.8148068) q[1];
sx q[1];
rz(-1.544373) q[1];
rz(-pi) q[2];
rz(0.88313734) q[3];
sx q[3];
rz(-2.4534907) q[3];
sx q[3];
rz(0.36206743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0110317) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(-1.4981184) q[2];
rz(2.2091852) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(1.9278661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907923) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(2.2687268) q[1];
sx q[1];
rz(-2.0582336) q[1];
sx q[1];
rz(2.2179926) q[1];
rz(-2.2070956) q[2];
sx q[2];
rz(-2.1515982) q[2];
sx q[2];
rz(1.2610255) q[2];
rz(-3.031293) q[3];
sx q[3];
rz(-0.58021373) q[3];
sx q[3];
rz(1.8662682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
