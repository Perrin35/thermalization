OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(4.1133147) q[1];
sx q[1];
rz(9.9749554) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101947) q[0];
sx q[0];
rz(-2.9203127) q[0];
sx q[0];
rz(1.6943323) q[0];
x q[1];
rz(-2.2416441) q[2];
sx q[2];
rz(-2.4953825) q[2];
sx q[2];
rz(-2.3044569) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3783592) q[1];
sx q[1];
rz(-2.2765056) q[1];
sx q[1];
rz(2.6325429) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4518277) q[3];
sx q[3];
rz(-2.527378) q[3];
sx q[3];
rz(0.14106942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(-2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(-2.7602592) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752276) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(-2.8711328) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(1.6289904) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0727901) q[0];
sx q[0];
rz(-1.8913942) q[0];
sx q[0];
rz(-0.94706236) q[0];
rz(1.5741882) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(1.1018745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7035547) q[1];
sx q[1];
rz(-2.2433271) q[1];
sx q[1];
rz(-0.1708252) q[1];
rz(-pi) q[2];
rz(2.9019722) q[3];
sx q[3];
rz(-1.8488956) q[3];
sx q[3];
rz(0.77277377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(2.846068) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(2.3957516) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5999818) q[0];
sx q[0];
rz(-2.9020502) q[0];
sx q[0];
rz(1.177686) q[0];
x q[1];
rz(-2.8638641) q[2];
sx q[2];
rz(-2.5979497) q[2];
sx q[2];
rz(-0.03253983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.44887603) q[1];
sx q[1];
rz(-1.8489031) q[1];
sx q[1];
rz(2.6443308) q[1];
rz(-1.8196705) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(2.4528465) q[2];
rz(3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.922309) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(-3.0396089) q[0];
rz(3.02137) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(2.8682958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7057987) q[0];
sx q[0];
rz(-1.4819643) q[0];
sx q[0];
rz(1.2466794) q[0];
rz(-pi) q[1];
rz(-0.69584537) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(2.3455182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4218688) q[1];
sx q[1];
rz(-0.27742741) q[1];
sx q[1];
rz(-1.006702) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9076505) q[3];
sx q[3];
rz(-0.66286874) q[3];
sx q[3];
rz(0.58594221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3499202) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824317) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(-3.058847) q[0];
rz(-0.67963183) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(-2.1544429) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521096) q[0];
sx q[0];
rz(-1.9728567) q[0];
sx q[0];
rz(-2.2228918) q[0];
x q[1];
rz(2.4977788) q[2];
sx q[2];
rz(-1.5632731) q[2];
sx q[2];
rz(1.8260337) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.134569) q[1];
sx q[1];
rz(-2.1937074) q[1];
sx q[1];
rz(-1.182215) q[1];
rz(-pi) q[2];
rz(1.1779285) q[3];
sx q[3];
rz(-1.6097027) q[3];
sx q[3];
rz(-1.9143357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(2.9303072) q[2];
rz(2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.485065) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(-0.791839) q[0];
rz(-0.99545288) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50514454) q[0];
sx q[0];
rz(-2.8594058) q[0];
sx q[0];
rz(1.6237153) q[0];
x q[1];
rz(-0.40462599) q[2];
sx q[2];
rz(-0.64794108) q[2];
sx q[2];
rz(2.4857156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4931603) q[1];
sx q[1];
rz(-1.5409768) q[1];
sx q[1];
rz(-1.0489419) q[1];
x q[2];
rz(-0.91026129) q[3];
sx q[3];
rz(-0.22781867) q[3];
sx q[3];
rz(-1.6662625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(0.77077579) q[2];
rz(-1.6714913) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8188748) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(-0.21559134) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(-0.0035704426) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63821793) q[0];
sx q[0];
rz(-2.5711381) q[0];
sx q[0];
rz(2.4128782) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9622757) q[2];
sx q[2];
rz(-2.2841095) q[2];
sx q[2];
rz(-1.6187402) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53686245) q[1];
sx q[1];
rz(-0.68568789) q[1];
sx q[1];
rz(-2.0525949) q[1];
rz(1.9332063) q[3];
sx q[3];
rz(-2.9299195) q[3];
sx q[3];
rz(3.1162804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-2.2214831) q[2];
rz(-2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6253117) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(0.40813804) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-0.67869854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40020254) q[0];
sx q[0];
rz(-0.42660248) q[0];
sx q[0];
rz(0.68047561) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0661725) q[2];
sx q[2];
rz(-2.3373211) q[2];
sx q[2];
rz(1.7240766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39290909) q[1];
sx q[1];
rz(-2.6137685) q[1];
sx q[1];
rz(1.4455568) q[1];
rz(-pi) q[2];
x q[2];
rz(2.802556) q[3];
sx q[3];
rz(-2.6936274) q[3];
sx q[3];
rz(2.8807004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(-1.142189) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1552102) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(2.881799) q[0];
rz(-2.4329176) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(2.6146467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5567112) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(2.7816539) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2692659) q[2];
sx q[2];
rz(-0.93062799) q[2];
sx q[2];
rz(0.38538853) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.770551) q[1];
sx q[1];
rz(-2.7275804) q[1];
sx q[1];
rz(0.71055926) q[1];
x q[2];
rz(-2.5832289) q[3];
sx q[3];
rz(-0.51279587) q[3];
sx q[3];
rz(0.59190291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-0.33716831) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-0.98544425) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-2.7808166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0246968) q[0];
sx q[0];
rz(-1.7381867) q[0];
sx q[0];
rz(-3.1162529) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.785727) q[2];
sx q[2];
rz(-1.4533918) q[2];
sx q[2];
rz(1.4431151) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39721397) q[1];
sx q[1];
rz(-0.75512868) q[1];
sx q[1];
rz(-0.20882864) q[1];
rz(1.4233227) q[3];
sx q[3];
rz(-1.7383988) q[3];
sx q[3];
rz(2.9180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(-1.7452314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5671134) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(1.7760989) q[2];
sx q[2];
rz(-1.4751954) q[2];
sx q[2];
rz(1.8196646) q[2];
rz(-0.75136649) q[3];
sx q[3];
rz(-2.1204227) q[3];
sx q[3];
rz(-1.942996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
