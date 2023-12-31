OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(6.0072748) q[0];
sx q[0];
rz(10.732565) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(-2.2059031) q[1];
sx q[1];
rz(1.5712665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21068621) q[0];
sx q[0];
rz(-1.9084198) q[0];
sx q[0];
rz(0.36436413) q[0];
x q[1];
rz(2.4351032) q[2];
sx q[2];
rz(-2.2405365) q[2];
sx q[2];
rz(-1.1342088) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1085514) q[1];
sx q[1];
rz(-1.3290977) q[1];
sx q[1];
rz(-0.34717314) q[1];
rz(2.0516112) q[3];
sx q[3];
rz(-2.7365723) q[3];
sx q[3];
rz(0.68457505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87542614) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(1.1323294) q[2];
rz(-1.4663565) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(2.9557513) q[0];
rz(-2.5813685) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-2.9247608) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8502055) q[0];
sx q[0];
rz(-0.7190401) q[0];
sx q[0];
rz(-2.015381) q[0];
rz(-2.2611513) q[2];
sx q[2];
rz(-2.464622) q[2];
sx q[2];
rz(1.134269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.617802) q[1];
sx q[1];
rz(-2.4317867) q[1];
sx q[1];
rz(1.8989423) q[1];
x q[2];
rz(0.87644491) q[3];
sx q[3];
rz(-2.090824) q[3];
sx q[3];
rz(-2.2976573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(-1.2878093) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(2.8365703) q[3];
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
rz(-1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(2.537354) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(-0.93260971) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9661449) q[0];
sx q[0];
rz(-1.6245337) q[0];
sx q[0];
rz(1.2530112) q[0];
rz(2.9595397) q[2];
sx q[2];
rz(-1.7943873) q[2];
sx q[2];
rz(1.455866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8285117) q[1];
sx q[1];
rz(-2.5701437) q[1];
sx q[1];
rz(-2.9213195) q[1];
rz(-pi) q[2];
rz(-2.1266537) q[3];
sx q[3];
rz(-1.1851289) q[3];
sx q[3];
rz(-2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(0.5422194) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820691) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(-0.50022593) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.4979699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1456137) q[0];
sx q[0];
rz(-1.0756452) q[0];
sx q[0];
rz(-0.33546319) q[0];
rz(1.285032) q[2];
sx q[2];
rz(-2.9432202) q[2];
sx q[2];
rz(-2.6464268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8343463) q[1];
sx q[1];
rz(-1.3386968) q[1];
sx q[1];
rz(2.8744065) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3718932) q[3];
sx q[3];
rz(-0.61004988) q[3];
sx q[3];
rz(-1.1806012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74636373) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(2.4397819) q[2];
rz(-2.3102405) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410626) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(-2.3262614) q[0];
rz(-1.5218081) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(2.0933847) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3807555) q[0];
sx q[0];
rz(-1.3073982) q[0];
sx q[0];
rz(-0.25816985) q[0];
rz(2.4900715) q[2];
sx q[2];
rz(-1.5610352) q[2];
sx q[2];
rz(-2.5196645) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8480307) q[1];
sx q[1];
rz(-1.3333496) q[1];
sx q[1];
rz(0.32229396) q[1];
rz(-pi) q[2];
rz(-1.0190373) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(0.90941959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(-1.099951) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0734171) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-0.90240479) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(3.0117603) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3467305) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(-2.5559588) q[0];
rz(-pi) q[1];
rz(-2.7733299) q[2];
sx q[2];
rz(-2.1149181) q[2];
sx q[2];
rz(0.95168176) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2643913) q[1];
sx q[1];
rz(-2.1346722) q[1];
sx q[1];
rz(2.0045723) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3744266) q[3];
sx q[3];
rz(-2.11103) q[3];
sx q[3];
rz(1.9062717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3123902) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(-0.20425805) q[2];
rz(1.9355109) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(-2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181353) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.6947421) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(2.4553305) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8780554) q[0];
sx q[0];
rz(-1.1362846) q[0];
sx q[0];
rz(2.6146019) q[0];
rz(-2.5574066) q[2];
sx q[2];
rz(-2.2224732) q[2];
sx q[2];
rz(2.3551031) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9069179) q[1];
sx q[1];
rz(-1.7559116) q[1];
sx q[1];
rz(0.075637416) q[1];
x q[2];
rz(0.42240123) q[3];
sx q[3];
rz(-1.8654612) q[3];
sx q[3];
rz(2.9871123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4454322) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(3.1398204) q[2];
rz(0.56162515) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(-1.6368438) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6034265) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(0.7810477) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.6400281) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9672464) q[0];
sx q[0];
rz(-1.5902728) q[0];
sx q[0];
rz(2.0635701) q[0];
x q[1];
rz(3.115032) q[2];
sx q[2];
rz(-1.6325258) q[2];
sx q[2];
rz(-2.314687) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0280684) q[1];
sx q[1];
rz(-1.6599732) q[1];
sx q[1];
rz(2.813617) q[1];
x q[2];
rz(-0.6286962) q[3];
sx q[3];
rz(-1.0843715) q[3];
sx q[3];
rz(1.4592255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7897196) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(2.7899172) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29031819) q[0];
sx q[0];
rz(-1.1997249) q[0];
sx q[0];
rz(-2.0593658) q[0];
rz(-2.7792764) q[2];
sx q[2];
rz(-3*pi/13) q[2];
sx q[2];
rz(-0.16513261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55191509) q[1];
sx q[1];
rz(-0.81112408) q[1];
sx q[1];
rz(-0.07304904) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44378186) q[3];
sx q[3];
rz(-3.0203331) q[3];
sx q[3];
rz(-1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7982771) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(0.19432755) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(2.1077572) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4170096) q[0];
sx q[0];
rz(-1.4405182) q[0];
sx q[0];
rz(-0.91086046) q[0];
x q[1];
rz(-2.6331484) q[2];
sx q[2];
rz(-2.2061081) q[2];
sx q[2];
rz(0.19257643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.918805) q[1];
sx q[1];
rz(-2.5606887) q[1];
sx q[1];
rz(-1.3851628) q[1];
rz(-2.8909573) q[3];
sx q[3];
rz(-2.4126629) q[3];
sx q[3];
rz(0.38754101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0620492) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-0.6357843) q[2];
rz(-2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.7059965) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(-0.031899115) q[2];
sx q[2];
rz(-2.1733641) q[2];
sx q[2];
rz(2.7410438) q[2];
rz(1.3676436) q[3];
sx q[3];
rz(-2.805134) q[3];
sx q[3];
rz(0.72611879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
