OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5306659) q[0];
sx q[0];
rz(-2.2025684) q[0];
sx q[0];
rz(0.0052069081) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(-1.985328) q[1];
sx q[1];
rz(-1.1896689) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.386728) q[0];
sx q[0];
rz(-2.0804188) q[0];
sx q[0];
rz(-1.1763563) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4283242) q[2];
sx q[2];
rz(-1.7490897) q[2];
sx q[2];
rz(1.1536191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8813821) q[1];
sx q[1];
rz(-1.8799026) q[1];
sx q[1];
rz(2.9865772) q[1];
rz(0.40640229) q[3];
sx q[3];
rz(-2.1592525) q[3];
sx q[3];
rz(-1.7455846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4253915) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7146724) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(-2.0479653) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(0.11322583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118606) q[0];
sx q[0];
rz(-1.9959873) q[0];
sx q[0];
rz(-3.0990764) q[0];
rz(-3.0683238) q[2];
sx q[2];
rz(-1.5511302) q[2];
sx q[2];
rz(-2.427223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6254239) q[1];
sx q[1];
rz(-0.58991573) q[1];
sx q[1];
rz(0.77074681) q[1];
rz(-pi) q[2];
rz(1.232997) q[3];
sx q[3];
rz(-0.63159734) q[3];
sx q[3];
rz(-2.256957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(-0.22228995) q[2];
rz(2.9120581) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(-0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753321) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(2.0942028) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(-3.0139794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532967) q[0];
sx q[0];
rz(-1.3971546) q[0];
sx q[0];
rz(2.2542473) q[0];
rz(-pi) q[1];
rz(2.2682842) q[2];
sx q[2];
rz(-1.0190939) q[2];
sx q[2];
rz(-1.973195) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0652005) q[1];
sx q[1];
rz(-0.25307357) q[1];
sx q[1];
rz(1.9051001) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94362887) q[3];
sx q[3];
rz(-1.3645932) q[3];
sx q[3];
rz(2.9507153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3601274) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(-2.8313417) q[2];
rz(-2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(-1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67868245) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(-0.15790766) q[0];
rz(2.7754916) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(-0.37240949) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5126257) q[0];
sx q[0];
rz(-0.86932875) q[0];
sx q[0];
rz(2.9823098) q[0];
rz(2.0580975) q[2];
sx q[2];
rz(-2.5022025) q[2];
sx q[2];
rz(2.1208178) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.403703) q[1];
sx q[1];
rz(-0.40778128) q[1];
sx q[1];
rz(-1.9562734) q[1];
rz(0.59471547) q[3];
sx q[3];
rz(-2.3200431) q[3];
sx q[3];
rz(-2.6298414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(-0.27238971) q[2];
rz(0.90879905) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.882778) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(2.0102665) q[0];
rz(-2.5436026) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(-2.4647443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15166053) q[0];
sx q[0];
rz(-2.6940072) q[0];
sx q[0];
rz(2.9582892) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9025181) q[2];
sx q[2];
rz(-0.74547807) q[2];
sx q[2];
rz(2.4649232) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.685806) q[1];
sx q[1];
rz(-1.5450486) q[1];
sx q[1];
rz(-2.8867433) q[1];
rz(2.7480514) q[3];
sx q[3];
rz(-1.593309) q[3];
sx q[3];
rz(-2.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(0.14906135) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1260219) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(1.4591249) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(-0.93313342) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3286896) q[0];
sx q[0];
rz(-2.0314386) q[0];
sx q[0];
rz(-2.8115389) q[0];
rz(1.2007347) q[2];
sx q[2];
rz(-1.0002631) q[2];
sx q[2];
rz(-0.42966118) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16490368) q[1];
sx q[1];
rz(-2.7737244) q[1];
sx q[1];
rz(-1.5771754) q[1];
rz(2.9643781) q[3];
sx q[3];
rz(-2.2209077) q[3];
sx q[3];
rz(-2.9110416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78559819) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(2.8549426) q[2];
rz(-2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048112415) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.6749143) q[0];
rz(1.02007) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(1.0010304) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0953656) q[0];
sx q[0];
rz(-1.2903178) q[0];
sx q[0];
rz(3.1026955) q[0];
rz(-0.83058968) q[2];
sx q[2];
rz(-2.6579654) q[2];
sx q[2];
rz(2.4110576) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.172799) q[1];
sx q[1];
rz(-1.4152923) q[1];
sx q[1];
rz(-1.96442) q[1];
x q[2];
rz(-0.27490297) q[3];
sx q[3];
rz(-2.4723845) q[3];
sx q[3];
rz(-1.6717403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.137407) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(3.0380847) q[2];
rz(0.59182709) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25933927) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(-0.6119588) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-2.7517095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.10605) q[0];
sx q[0];
rz(-2.4379726) q[0];
sx q[0];
rz(2.5662867) q[0];
rz(-pi) q[1];
rz(-2.6340918) q[2];
sx q[2];
rz(-1.9365053) q[2];
sx q[2];
rz(0.65069228) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0695614) q[1];
sx q[1];
rz(-1.8998635) q[1];
sx q[1];
rz(-2.0272994) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89160664) q[3];
sx q[3];
rz(-1.3636175) q[3];
sx q[3];
rz(0.87219119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(2.4273382) q[2];
rz(0.8977302) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(-0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(1.6927239) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(2.3642448) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(2.8651967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9696635) q[0];
sx q[0];
rz(-1.164325) q[0];
sx q[0];
rz(1.7448234) q[0];
rz(1.3771636) q[2];
sx q[2];
rz(-1.7948705) q[2];
sx q[2];
rz(-0.60806882) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8564057) q[1];
sx q[1];
rz(-1.7289843) q[1];
sx q[1];
rz(-1.3473947) q[1];
rz(-1.8180088) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(1.1288527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(0.12750553) q[2];
rz(0.036711983) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7413095) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(0.35274831) q[0];
rz(2.5648975) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(-1.030285) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.661783) q[0];
sx q[0];
rz(-2.3942238) q[0];
sx q[0];
rz(2.8519147) q[0];
rz(-pi) q[1];
rz(-2.7416517) q[2];
sx q[2];
rz(-2.415014) q[2];
sx q[2];
rz(0.71619294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.13223091) q[1];
sx q[1];
rz(-1.2538326) q[1];
sx q[1];
rz(-0.46225458) q[1];
rz(-0.45566166) q[3];
sx q[3];
rz(-0.66484287) q[3];
sx q[3];
rz(-1.957422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0673922) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-0.34118787) q[2];
rz(-2.4009005) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(3.0497131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13070233) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(-2.7080766) q[2];
sx q[2];
rz(-1.3939861) q[2];
sx q[2];
rz(1.3409333) q[2];
rz(-2.1445198) q[3];
sx q[3];
rz(-2.1763747) q[3];
sx q[3];
rz(2.5381266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
