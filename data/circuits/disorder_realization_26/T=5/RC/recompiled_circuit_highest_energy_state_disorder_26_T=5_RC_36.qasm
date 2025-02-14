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
rz(-2.969279) q[0];
sx q[0];
rz(-0.20792374) q[0];
sx q[0];
rz(-1.2907668) q[0];
rz(-2.9895904) q[1];
sx q[1];
rz(-0.64270371) q[1];
sx q[1];
rz(-0.42833498) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5458459) q[0];
sx q[0];
rz(-1.6753249) q[0];
sx q[0];
rz(1.680611) q[0];
rz(-pi) q[1];
x q[1];
rz(1.27698) q[2];
sx q[2];
rz(-0.91962469) q[2];
sx q[2];
rz(3.0319954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1312771) q[1];
sx q[1];
rz(-0.73956623) q[1];
sx q[1];
rz(-0.063994813) q[1];
x q[2];
rz(-2.9228404) q[3];
sx q[3];
rz(-1.7667337) q[3];
sx q[3];
rz(2.4234555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20994818) q[2];
sx q[2];
rz(-2.1817709) q[2];
sx q[2];
rz(1.2098562) q[2];
rz(-0.079553902) q[3];
sx q[3];
rz(-0.39593655) q[3];
sx q[3];
rz(2.615926) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88456589) q[0];
sx q[0];
rz(-2.631812) q[0];
sx q[0];
rz(2.7575745) q[0];
rz(-2.2513023) q[1];
sx q[1];
rz(-1.236981) q[1];
sx q[1];
rz(0.30337897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6174406) q[0];
sx q[0];
rz(-2.5684803) q[0];
sx q[0];
rz(-1.9277431) q[0];
x q[1];
rz(1.0531002) q[2];
sx q[2];
rz(-1.7974904) q[2];
sx q[2];
rz(2.0402619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5785256) q[1];
sx q[1];
rz(-1.6614698) q[1];
sx q[1];
rz(3.0791326) q[1];
x q[2];
rz(2.5796765) q[3];
sx q[3];
rz(-0.73560916) q[3];
sx q[3];
rz(0.058678415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38637912) q[2];
sx q[2];
rz(-1.1073802) q[2];
sx q[2];
rz(-2.4066822) q[2];
rz(0.84878659) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(0.36062226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3291149) q[0];
sx q[0];
rz(-0.37549967) q[0];
sx q[0];
rz(3.062881) q[0];
rz(-2.3178237) q[1];
sx q[1];
rz(-2.0562101) q[1];
sx q[1];
rz(1.8471921) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9322008) q[0];
sx q[0];
rz(-3.0552312) q[0];
sx q[0];
rz(-2.3461653) q[0];
x q[1];
rz(-2.6316775) q[2];
sx q[2];
rz(-1.0140061) q[2];
sx q[2];
rz(2.315588) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8405806) q[1];
sx q[1];
rz(-2.8361179) q[1];
sx q[1];
rz(-2.7423647) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23842509) q[3];
sx q[3];
rz(-0.88078558) q[3];
sx q[3];
rz(-2.391614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89982975) q[2];
sx q[2];
rz(-2.7624942) q[2];
sx q[2];
rz(2.3251593) q[2];
rz(-1.624931) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(0.70037705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65275943) q[0];
sx q[0];
rz(-2.3214898) q[0];
sx q[0];
rz(-2.5732727) q[0];
rz(-1.142451) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(-1.6894587) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84684337) q[0];
sx q[0];
rz(-0.94943014) q[0];
sx q[0];
rz(-1.2970379) q[0];
rz(1.6644888) q[2];
sx q[2];
rz(-1.7498651) q[2];
sx q[2];
rz(-0.25824091) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8206222) q[1];
sx q[1];
rz(-1.2577211) q[1];
sx q[1];
rz(-2.2105107) q[1];
rz(0.178517) q[3];
sx q[3];
rz(-1.2367289) q[3];
sx q[3];
rz(0.26135437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51231724) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(-3.0647035) q[2];
rz(2.7422089) q[3];
sx q[3];
rz(-0.9136343) q[3];
sx q[3];
rz(-0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9870616) q[0];
sx q[0];
rz(-2.7847544) q[0];
sx q[0];
rz(-2.8416908) q[0];
rz(-1.9918282) q[1];
sx q[1];
rz(-1.0118142) q[1];
sx q[1];
rz(2.6531175) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73672527) q[0];
sx q[0];
rz(-2.953385) q[0];
sx q[0];
rz(1.7467278) q[0];
rz(-pi) q[1];
rz(-0.68814338) q[2];
sx q[2];
rz(-2.4240085) q[2];
sx q[2];
rz(-0.46844278) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3255182) q[1];
sx q[1];
rz(-1.5897004) q[1];
sx q[1];
rz(-0.58362785) q[1];
rz(-pi) q[2];
rz(-2.4998922) q[3];
sx q[3];
rz(-1.7880758) q[3];
sx q[3];
rz(-3.1171796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7601295) q[2];
sx q[2];
rz(-3.0035786) q[2];
sx q[2];
rz(1.4786973) q[2];
rz(-2.0315157) q[3];
sx q[3];
rz(-2.1434982) q[3];
sx q[3];
rz(-2.4983675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396624) q[0];
sx q[0];
rz(-1.9597541) q[0];
sx q[0];
rz(2.8231743) q[0];
rz(3.1069801) q[1];
sx q[1];
rz(-2.5739659) q[1];
sx q[1];
rz(2.6908223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0552669) q[0];
sx q[0];
rz(-2.0811715) q[0];
sx q[0];
rz(-0.016252131) q[0];
x q[1];
rz(-2.320596) q[2];
sx q[2];
rz(-0.85047532) q[2];
sx q[2];
rz(-0.79325097) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3035982) q[1];
sx q[1];
rz(-2.384409) q[1];
sx q[1];
rz(0.96214575) q[1];
x q[2];
rz(-0.3698753) q[3];
sx q[3];
rz(-2.2910017) q[3];
sx q[3];
rz(-2.6323619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0938809) q[2];
sx q[2];
rz(-0.1897976) q[2];
sx q[2];
rz(0.06165687) q[2];
rz(-0.098585248) q[3];
sx q[3];
rz(-0.74461377) q[3];
sx q[3];
rz(-0.70257598) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3848569) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(3.1016438) q[0];
rz(-0.7705676) q[1];
sx q[1];
rz(-0.71208411) q[1];
sx q[1];
rz(0.22824731) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6421445) q[0];
sx q[0];
rz(-0.11760437) q[0];
sx q[0];
rz(2.2884877) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1399916) q[2];
sx q[2];
rz(-1.8769771) q[2];
sx q[2];
rz(-1.209335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3769987) q[1];
sx q[1];
rz(-1.5176763) q[1];
sx q[1];
rz(-0.082100987) q[1];
x q[2];
rz(0.41577783) q[3];
sx q[3];
rz(-2.1630187) q[3];
sx q[3];
rz(1.3982738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.063529) q[2];
sx q[2];
rz(-1.0588131) q[2];
sx q[2];
rz(0.25827363) q[2];
rz(-0.20049788) q[3];
sx q[3];
rz(-2.343488) q[3];
sx q[3];
rz(0.18331461) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9789326) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(-0.3072511) q[0];
rz(-1.2112674) q[1];
sx q[1];
rz(-0.4937506) q[1];
sx q[1];
rz(0.58120751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8756105) q[0];
sx q[0];
rz(-1.6020244) q[0];
sx q[0];
rz(-1.9755548) q[0];
rz(-pi) q[1];
rz(2.4318784) q[2];
sx q[2];
rz(-1.6953371) q[2];
sx q[2];
rz(-2.1309851) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30347428) q[1];
sx q[1];
rz(-1.3708769) q[1];
sx q[1];
rz(-2.3568627) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56116207) q[3];
sx q[3];
rz(-0.61120874) q[3];
sx q[3];
rz(1.2964013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6792949) q[2];
sx q[2];
rz(-2.0079948) q[2];
sx q[2];
rz(3.1104258) q[2];
rz(0.22552414) q[3];
sx q[3];
rz(-1.3616819) q[3];
sx q[3];
rz(-2.1112736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0533326) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(2.4326676) q[0];
rz(0.21253474) q[1];
sx q[1];
rz(-2.1268851) q[1];
sx q[1];
rz(-0.85740352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5564726) q[0];
sx q[0];
rz(-1.5758638) q[0];
sx q[0];
rz(1.7005672) q[0];
x q[1];
rz(-1.2456263) q[2];
sx q[2];
rz(-2.0813) q[2];
sx q[2];
rz(2.4028962) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6826806) q[1];
sx q[1];
rz(-2.4574728) q[1];
sx q[1];
rz(0.64351179) q[1];
rz(-pi) q[2];
rz(0.25257381) q[3];
sx q[3];
rz(-0.71306224) q[3];
sx q[3];
rz(-2.3424462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4400441) q[2];
sx q[2];
rz(-2.7893119) q[2];
sx q[2];
rz(2.3360543) q[2];
rz(-2.7696179) q[3];
sx q[3];
rz(-1.493908) q[3];
sx q[3];
rz(2.7443547) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1086248) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(0.97880542) q[0];
rz(-2.4188614) q[1];
sx q[1];
rz(-1.1284072) q[1];
sx q[1];
rz(-0.61789787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5586097) q[0];
sx q[0];
rz(-1.2961868) q[0];
sx q[0];
rz(-2.3586169) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63510908) q[2];
sx q[2];
rz(-2.375787) q[2];
sx q[2];
rz(2.2857411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5717575) q[1];
sx q[1];
rz(-1.7392842) q[1];
sx q[1];
rz(1.38901) q[1];
rz(0.083656351) q[3];
sx q[3];
rz(-1.6550875) q[3];
sx q[3];
rz(2.313446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2021947) q[2];
sx q[2];
rz(-1.3556182) q[2];
sx q[2];
rz(-0.00016577684) q[2];
rz(2.557142) q[3];
sx q[3];
rz(-2.1355459) q[3];
sx q[3];
rz(2.5463026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38681876) q[0];
sx q[0];
rz(-1.5712354) q[0];
sx q[0];
rz(-1.5729217) q[0];
rz(-1.3407002) q[1];
sx q[1];
rz(-1.1005713) q[1];
sx q[1];
rz(1.5250199) q[1];
rz(-0.93201119) q[2];
sx q[2];
rz(-1.1975653) q[2];
sx q[2];
rz(-1.4690659) q[2];
rz(0.67047337) q[3];
sx q[3];
rz(-2.5994876) q[3];
sx q[3];
rz(2.2710298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
