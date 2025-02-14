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
rz(0.33557284) q[0];
sx q[0];
rz(-1.4541452) q[0];
sx q[0];
rz(-2.1079221) q[0];
rz(1.4312862) q[1];
sx q[1];
rz(-0.55748504) q[1];
sx q[1];
rz(2.1318336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97573853) q[0];
sx q[0];
rz(-1.3814808) q[0];
sx q[0];
rz(-3.0136787) q[0];
rz(-1.4883243) q[2];
sx q[2];
rz(-1.5745275) q[2];
sx q[2];
rz(-0.53860215) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7518471) q[1];
sx q[1];
rz(-1.455014) q[1];
sx q[1];
rz(0.16645653) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5752154) q[3];
sx q[3];
rz(-1.4285681) q[3];
sx q[3];
rz(2.9517236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3681763) q[2];
sx q[2];
rz(-0.877031) q[2];
sx q[2];
rz(-2.3167493) q[2];
rz(-0.33251897) q[3];
sx q[3];
rz(-0.85086346) q[3];
sx q[3];
rz(2.6578145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2291718) q[0];
sx q[0];
rz(-3.0569172) q[0];
sx q[0];
rz(0.53572768) q[0];
rz(2.3748705) q[1];
sx q[1];
rz(-1.7886536) q[1];
sx q[1];
rz(1.7492693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77477876) q[0];
sx q[0];
rz(-2.5486425) q[0];
sx q[0];
rz(-2.2137027) q[0];
x q[1];
rz(0.34113348) q[2];
sx q[2];
rz(-1.5798414) q[2];
sx q[2];
rz(1.1760528) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0619095) q[1];
sx q[1];
rz(-2.3259729) q[1];
sx q[1];
rz(0.029600005) q[1];
rz(-pi) q[2];
rz(1.6757319) q[3];
sx q[3];
rz(-1.4426878) q[3];
sx q[3];
rz(-0.48612938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2367737) q[2];
sx q[2];
rz(-2.3084013) q[2];
sx q[2];
rz(1.0986249) q[2];
rz(-0.83313471) q[3];
sx q[3];
rz(-1.5435217) q[3];
sx q[3];
rz(1.0540849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5611834) q[0];
sx q[0];
rz(-1.9920749) q[0];
sx q[0];
rz(-2.4533601) q[0];
rz(1.2511823) q[1];
sx q[1];
rz(-2.4871608) q[1];
sx q[1];
rz(1.9293264) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1935342) q[0];
sx q[0];
rz(-2.2431264) q[0];
sx q[0];
rz(2.1904519) q[0];
x q[1];
rz(-2.9040292) q[2];
sx q[2];
rz(-1.8389987) q[2];
sx q[2];
rz(-2.149127) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1656629) q[1];
sx q[1];
rz(-1.5651795) q[1];
sx q[1];
rz(-0.65896874) q[1];
rz(1.4914277) q[3];
sx q[3];
rz(-1.5659596) q[3];
sx q[3];
rz(1.2629379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4091829) q[2];
sx q[2];
rz(-2.1467291) q[2];
sx q[2];
rz(-2.691972) q[2];
rz(0.85363394) q[3];
sx q[3];
rz(-0.64695224) q[3];
sx q[3];
rz(-2.4322815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7483826) q[0];
sx q[0];
rz(-2.125183) q[0];
sx q[0];
rz(2.7384695) q[0];
rz(-0.77278167) q[1];
sx q[1];
rz(-1.2330331) q[1];
sx q[1];
rz(0.82537878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6306303) q[0];
sx q[0];
rz(-1.2494506) q[0];
sx q[0];
rz(2.3031631) q[0];
rz(2.4140777) q[2];
sx q[2];
rz(-3.10559) q[2];
sx q[2];
rz(1.9787479) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.17887184) q[1];
sx q[1];
rz(-0.58607465) q[1];
sx q[1];
rz(1.7078215) q[1];
rz(-pi) q[2];
rz(-2.585936) q[3];
sx q[3];
rz(-0.98816493) q[3];
sx q[3];
rz(-0.61911303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0394502) q[2];
sx q[2];
rz(-1.8723698) q[2];
sx q[2];
rz(0.5086745) q[2];
rz(1.2447119) q[3];
sx q[3];
rz(-1.0735984) q[3];
sx q[3];
rz(-1.3332453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73388571) q[0];
sx q[0];
rz(-1.5835967) q[0];
sx q[0];
rz(0.35633126) q[0];
rz(-1.7286667) q[1];
sx q[1];
rz(-1.3099542) q[1];
sx q[1];
rz(-1.2948571) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6083487) q[0];
sx q[0];
rz(-2.6496467) q[0];
sx q[0];
rz(-2.66376) q[0];
rz(-pi) q[1];
rz(-1.4003428) q[2];
sx q[2];
rz(-0.91636412) q[2];
sx q[2];
rz(2.7026388) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7562189) q[1];
sx q[1];
rz(-1.1750147) q[1];
sx q[1];
rz(-1.8204921) q[1];
rz(-pi) q[2];
rz(1.0255085) q[3];
sx q[3];
rz(-1.9615615) q[3];
sx q[3];
rz(-2.7440967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20092189) q[2];
sx q[2];
rz(-1.7565497) q[2];
sx q[2];
rz(-2.5015639) q[2];
rz(2.707543) q[3];
sx q[3];
rz(-2.1899352) q[3];
sx q[3];
rz(2.0210338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5459179) q[0];
sx q[0];
rz(-0.43788236) q[0];
sx q[0];
rz(-0.75089279) q[0];
rz(-0.76260507) q[1];
sx q[1];
rz(-1.4354939) q[1];
sx q[1];
rz(-2.6079752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55706829) q[0];
sx q[0];
rz(-1.8584588) q[0];
sx q[0];
rz(-2.8329222) q[0];
rz(-pi) q[1];
rz(-2.3377588) q[2];
sx q[2];
rz(-1.2495923) q[2];
sx q[2];
rz(-0.28576947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68249629) q[1];
sx q[1];
rz(-0.67855103) q[1];
sx q[1];
rz(1.7595923) q[1];
rz(-pi) q[2];
rz(-0.814416) q[3];
sx q[3];
rz(-2.4049063) q[3];
sx q[3];
rz(1.8131922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0642455) q[2];
sx q[2];
rz(-2.7392445) q[2];
sx q[2];
rz(-2.2557491) q[2];
rz(-1.3028076) q[3];
sx q[3];
rz(-1.9582483) q[3];
sx q[3];
rz(-2.6634789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.9293905) q[0];
sx q[0];
rz(-2.2181856) q[0];
sx q[0];
rz(2.5481664) q[0];
rz(1.5133096) q[1];
sx q[1];
rz(-1.1791041) q[1];
sx q[1];
rz(-2.5679307) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5705171) q[0];
sx q[0];
rz(-0.62643753) q[0];
sx q[0];
rz(-2.4721739) q[0];
rz(-pi) q[1];
rz(-1.8183579) q[2];
sx q[2];
rz(-1.0070966) q[2];
sx q[2];
rz(-0.42681387) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.03790126) q[1];
sx q[1];
rz(-1.173844) q[1];
sx q[1];
rz(1.3034225) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6214293) q[3];
sx q[3];
rz(-1.4922499) q[3];
sx q[3];
rz(2.1861064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4297428) q[2];
sx q[2];
rz(-1.9387127) q[2];
sx q[2];
rz(-1.708606) q[2];
rz(-1.1687219) q[3];
sx q[3];
rz(-0.43787268) q[3];
sx q[3];
rz(1.5168813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2699921) q[0];
sx q[0];
rz(-0.90684909) q[0];
sx q[0];
rz(2.9378743) q[0];
rz(1.3153971) q[1];
sx q[1];
rz(-1.306465) q[1];
sx q[1];
rz(-1.809583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4341605) q[0];
sx q[0];
rz(-0.95904175) q[0];
sx q[0];
rz(2.3513398) q[0];
rz(2.2384127) q[2];
sx q[2];
rz(-1.1051264) q[2];
sx q[2];
rz(0.47448928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9086919) q[1];
sx q[1];
rz(-0.76592839) q[1];
sx q[1];
rz(2.5604064) q[1];
x q[2];
rz(-2.2409954) q[3];
sx q[3];
rz(-0.20684563) q[3];
sx q[3];
rz(-1.8781917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77067644) q[2];
sx q[2];
rz(-3.0597661) q[2];
sx q[2];
rz(-1.3814242) q[2];
rz(-2.5904739) q[3];
sx q[3];
rz(-1.328732) q[3];
sx q[3];
rz(-0.74530017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119174) q[0];
sx q[0];
rz(-1.6918007) q[0];
sx q[0];
rz(-1.2592738) q[0];
rz(-1.7970386) q[1];
sx q[1];
rz(-2.5090736) q[1];
sx q[1];
rz(-1.9154027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38905242) q[0];
sx q[0];
rz(-1.2864094) q[0];
sx q[0];
rz(0.14099462) q[0];
x q[1];
rz(2.0001423) q[2];
sx q[2];
rz(-1.9252021) q[2];
sx q[2];
rz(0.62858519) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3873432) q[1];
sx q[1];
rz(-2.3918536) q[1];
sx q[1];
rz(3.0453338) q[1];
rz(-1.0223844) q[3];
sx q[3];
rz(-2.5520241) q[3];
sx q[3];
rz(1.1353253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5598477) q[2];
sx q[2];
rz(-2.2793844) q[2];
sx q[2];
rz(-1.5398514) q[2];
rz(-1.3567443) q[3];
sx q[3];
rz(-0.53250766) q[3];
sx q[3];
rz(-1.9073568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.27720472) q[0];
sx q[0];
rz(-1.5978403) q[0];
sx q[0];
rz(-3.0371015) q[0];
rz(1.744572) q[1];
sx q[1];
rz(-1.3518159) q[1];
sx q[1];
rz(1.6207961) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049690954) q[0];
sx q[0];
rz(-2.3713863) q[0];
sx q[0];
rz(0.30253221) q[0];
rz(-pi) q[1];
rz(2.7149523) q[2];
sx q[2];
rz(-2.7707515) q[2];
sx q[2];
rz(-1.1877354) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4222153) q[1];
sx q[1];
rz(-1.6746192) q[1];
sx q[1];
rz(0.28909282) q[1];
rz(-pi) q[2];
x q[2];
rz(0.057985882) q[3];
sx q[3];
rz(-0.99000217) q[3];
sx q[3];
rz(0.58422663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70667679) q[2];
sx q[2];
rz(-2.2997663) q[2];
sx q[2];
rz(-1.8527156) q[2];
rz(2.1364818) q[3];
sx q[3];
rz(-2.0821327) q[3];
sx q[3];
rz(-0.10678261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.91916753) q[0];
sx q[0];
rz(-1.5594302) q[0];
sx q[0];
rz(2.9710309) q[0];
rz(-2.2617321) q[1];
sx q[1];
rz(-2.6251371) q[1];
sx q[1];
rz(0.62216204) q[1];
rz(1.5063029) q[2];
sx q[2];
rz(-1.0721102) q[2];
sx q[2];
rz(-0.67244844) q[2];
rz(-2.5414657) q[3];
sx q[3];
rz(-1.5059581) q[3];
sx q[3];
rz(-3.0653421) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
