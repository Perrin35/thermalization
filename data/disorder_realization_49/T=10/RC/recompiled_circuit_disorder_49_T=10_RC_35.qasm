OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61092678) q[0];
sx q[0];
rz(-0.93902421) q[0];
sx q[0];
rz(-0.0052069081) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(-1.985328) q[1];
sx q[1];
rz(-1.1896689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016276377) q[0];
sx q[0];
rz(-1.9128886) q[0];
sx q[0];
rz(-2.5972511) q[0];
x q[1];
rz(2.4742545) q[2];
sx q[2];
rz(-0.22775209) q[2];
sx q[2];
rz(1.3078794) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78566879) q[1];
sx q[1];
rz(-2.7969116) q[1];
sx q[1];
rz(1.1204526) q[1];
x q[2];
rz(-0.40640229) q[3];
sx q[3];
rz(-0.98234017) q[3];
sx q[3];
rz(-1.7455846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71620119) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(0.5973967) q[2];
rz(-1.776009) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7146724) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(-2.0479653) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(-0.11322583) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118606) q[0];
sx q[0];
rz(-1.9959873) q[0];
sx q[0];
rz(3.0990764) q[0];
x q[1];
rz(3.0683238) q[2];
sx q[2];
rz(-1.5904625) q[2];
sx q[2];
rz(0.71436963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3791703) q[1];
sx q[1];
rz(-1.9813073) q[1];
sx q[1];
rz(-1.1344086) q[1];
x q[2];
rz(0.96666386) q[3];
sx q[3];
rz(-1.767744) q[3];
sx q[3];
rz(0.40991022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(2.9193027) q[2];
rz(-2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753321) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(3.0139794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.815044) q[0];
sx q[0];
rz(-0.70171261) q[0];
sx q[0];
rz(1.2998507) q[0];
rz(-pi) q[1];
rz(-0.87330841) q[2];
sx q[2];
rz(-1.0190939) q[2];
sx q[2];
rz(-1.973195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0652005) q[1];
sx q[1];
rz(-0.25307357) q[1];
sx q[1];
rz(-1.9051001) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25281275) q[3];
sx q[3];
rz(-2.1826934) q[3];
sx q[3];
rz(1.52724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7814653) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(0.310251) q[2];
rz(0.34960738) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.7278956) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-0.15790766) q[0];
rz(-2.7754916) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(0.37240949) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87293738) q[0];
sx q[0];
rz(-0.71632179) q[0];
sx q[0];
rz(1.3852081) q[0];
rz(-pi) q[1];
rz(-1.0834951) q[2];
sx q[2];
rz(-2.5022025) q[2];
sx q[2];
rz(-1.0207748) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9875033) q[1];
sx q[1];
rz(-1.9470864) q[1];
sx q[1];
rz(-2.980568) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5468772) q[3];
sx q[3];
rz(-0.82154951) q[3];
sx q[3];
rz(0.51175129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(-2.8692029) q[2];
rz(-2.2327936) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(1.1313261) q[0];
rz(0.5979901) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(0.67684832) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2535431) q[0];
sx q[0];
rz(-1.6497668) q[0];
sx q[0];
rz(-2.700564) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23907451) q[2];
sx q[2];
rz(-2.3961146) q[2];
sx q[2];
rz(-2.4649232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21351335) q[1];
sx q[1];
rz(-2.885474) q[1];
sx q[1];
rz(-3.03979) q[1];
x q[2];
rz(-2.7480514) q[3];
sx q[3];
rz(-1.593309) q[3];
sx q[3];
rz(-0.81641203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1375492) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(-0.14906135) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(-2.1242274) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.4591249) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(2.2084592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81290302) q[0];
sx q[0];
rz(-2.0314386) q[0];
sx q[0];
rz(-0.33005379) q[0];
rz(-pi) q[1];
rz(-0.51322333) q[2];
sx q[2];
rz(-0.66868082) q[2];
sx q[2];
rz(-0.19323397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.976689) q[1];
sx q[1];
rz(-2.7737244) q[1];
sx q[1];
rz(1.5771754) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.798614) q[3];
sx q[3];
rz(-2.4711547) q[3];
sx q[3];
rz(-3.0844641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(0.28665001) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048112415) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.4666784) q[0];
rz(2.1215227) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(-2.1405623) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2350378) q[0];
sx q[0];
rz(-0.28309238) q[0];
sx q[0];
rz(-1.7049768) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34044388) q[2];
sx q[2];
rz(-1.2203487) q[2];
sx q[2];
rz(-3.0710789) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8965473) q[1];
sx q[1];
rz(-2.7198615) q[1];
sx q[1];
rz(1.1827724) q[1];
rz(-1.3592968) q[3];
sx q[3];
rz(-0.930951) q[3];
sx q[3];
rz(-1.1246455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.137407) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(0.10350791) q[2];
rz(2.5497656) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25933927) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(2.7517095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.10605) q[0];
sx q[0];
rz(-2.4379726) q[0];
sx q[0];
rz(-2.5662867) q[0];
x q[1];
rz(-1.9837603) q[2];
sx q[2];
rz(-1.0997084) q[2];
sx q[2];
rz(0.72380356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4854359) q[1];
sx q[1];
rz(-1.1404783) q[1];
sx q[1];
rz(0.3635316) q[1];
rz(-pi) q[2];
rz(0.89160664) q[3];
sx q[3];
rz(-1.3636175) q[3];
sx q[3];
rz(-0.87219119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-2.4273382) q[2];
rz(2.2438625) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6927239) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(2.3642448) q[0];
rz(2.2946987) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(2.8651967) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9696635) q[0];
sx q[0];
rz(-1.164325) q[0];
sx q[0];
rz(-1.3967692) q[0];
x q[1];
rz(-1.7644291) q[2];
sx q[2];
rz(-1.7948705) q[2];
sx q[2];
rz(2.5335238) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8917577) q[1];
sx q[1];
rz(-1.7913622) q[1];
sx q[1];
rz(-2.9794429) q[1];
rz(-2.7941462) q[3];
sx q[3];
rz(-0.63784079) q[3];
sx q[3];
rz(-0.70536246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85236621) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(-0.12750553) q[2];
rz(-0.036711983) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7413095) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(0.35274831) q[0];
rz(-0.57669512) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(2.1113077) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.83537) q[0];
sx q[0];
rz(-1.3754002) q[0];
sx q[0];
rz(0.72619254) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39994098) q[2];
sx q[2];
rz(-2.415014) q[2];
sx q[2];
rz(-0.71619294) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5926338) q[1];
sx q[1];
rz(-1.1332129) q[1];
sx q[1];
rz(1.9220819) q[1];
rz(-pi) q[2];
rz(0.45566166) q[3];
sx q[3];
rz(-0.66484287) q[3];
sx q[3];
rz(1.957422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(0.34118787) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13070233) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(-1.7651991) q[2];
sx q[2];
rz(-1.1444848) q[2];
sx q[2];
rz(-0.3111006) q[2];
rz(-0.99707281) q[3];
sx q[3];
rz(-0.96521796) q[3];
sx q[3];
rz(-0.60346606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
