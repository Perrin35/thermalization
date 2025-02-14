OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0134861) q[0];
sx q[0];
rz(-0.82733893) q[0];
sx q[0];
rz(-1.3878571) q[0];
rz(-2.2881621) q[1];
sx q[1];
rz(-2.7397459) q[1];
sx q[1];
rz(-2.0274577) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3281109) q[0];
sx q[0];
rz(-0.8719647) q[0];
sx q[0];
rz(2.8100886) q[0];
rz(-pi) q[1];
rz(-2.4365455) q[2];
sx q[2];
rz(-1.4234666) q[2];
sx q[2];
rz(1.6697869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88190313) q[1];
sx q[1];
rz(-1.5659386) q[1];
sx q[1];
rz(1.5905981) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57350343) q[3];
sx q[3];
rz(-2.5636004) q[3];
sx q[3];
rz(-1.0740394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9806597) q[2];
sx q[2];
rz(-0.42155835) q[2];
sx q[2];
rz(1.6932311) q[2];
rz(2.8990922) q[3];
sx q[3];
rz(-0.91553965) q[3];
sx q[3];
rz(-1.2601132) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5777957) q[0];
sx q[0];
rz(-1.8255434) q[0];
sx q[0];
rz(1.0003566) q[0];
rz(0.73161221) q[1];
sx q[1];
rz(-1.4294521) q[1];
sx q[1];
rz(-1.9614722) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6574616) q[0];
sx q[0];
rz(-3.0041782) q[0];
sx q[0];
rz(0.71533393) q[0];
rz(-pi) q[1];
rz(3.1077049) q[2];
sx q[2];
rz(-1.8485957) q[2];
sx q[2];
rz(-0.44521618) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48881179) q[1];
sx q[1];
rz(-3.0074357) q[1];
sx q[1];
rz(1.556788) q[1];
x q[2];
rz(1.8537221) q[3];
sx q[3];
rz(-1.3045505) q[3];
sx q[3];
rz(0.11582813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0483027) q[2];
sx q[2];
rz(-0.90128171) q[2];
sx q[2];
rz(-1.0789336) q[2];
rz(-0.087873936) q[3];
sx q[3];
rz(-1.7885957) q[3];
sx q[3];
rz(2.8972076) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3946149) q[0];
sx q[0];
rz(-1.1629539) q[0];
sx q[0];
rz(1.278247) q[0];
rz(2.6784015) q[1];
sx q[1];
rz(-1.3563503) q[1];
sx q[1];
rz(0.82122222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6919977) q[0];
sx q[0];
rz(-2.2173081) q[0];
sx q[0];
rz(-2.9086935) q[0];
rz(0.37126705) q[2];
sx q[2];
rz(-1.5131931) q[2];
sx q[2];
rz(2.0774373) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48135172) q[1];
sx q[1];
rz(-1.4235745) q[1];
sx q[1];
rz(-0.46105701) q[1];
x q[2];
rz(2.3825112) q[3];
sx q[3];
rz(-1.0831175) q[3];
sx q[3];
rz(1.9239192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.897573) q[2];
sx q[2];
rz(-1.5622746) q[2];
sx q[2];
rz(2.139034) q[2];
rz(1.9133441) q[3];
sx q[3];
rz(-1.4017665) q[3];
sx q[3];
rz(1.4209411) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9035852) q[0];
sx q[0];
rz(-2.0915732) q[0];
sx q[0];
rz(2.8557657) q[0];
rz(-2.8803275) q[1];
sx q[1];
rz(-1.8087872) q[1];
sx q[1];
rz(0.030698311) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8809533) q[0];
sx q[0];
rz(-0.27361037) q[0];
sx q[0];
rz(0.27283313) q[0];
rz(-pi) q[1];
rz(0.89392406) q[2];
sx q[2];
rz(-1.5090164) q[2];
sx q[2];
rz(0.0047574818) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0568119) q[1];
sx q[1];
rz(-1.1263945) q[1];
sx q[1];
rz(0.2803352) q[1];
rz(1.8844344) q[3];
sx q[3];
rz(-1.6480069) q[3];
sx q[3];
rz(-1.9208637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3280481) q[2];
sx q[2];
rz(-1.8009461) q[2];
sx q[2];
rz(-3.1385341) q[2];
rz(2.3004153) q[3];
sx q[3];
rz(-2.5523461) q[3];
sx q[3];
rz(0.35471788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-2.4641814) q[0];
sx q[0];
rz(-2.2941636) q[0];
sx q[0];
rz(-1.4671951) q[0];
rz(-1.6293619) q[1];
sx q[1];
rz(-0.96577516) q[1];
sx q[1];
rz(-0.95219749) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.061082) q[0];
sx q[0];
rz(-3.0709776) q[0];
sx q[0];
rz(0.10297601) q[0];
rz(-pi) q[1];
rz(0.69779863) q[2];
sx q[2];
rz(-1.2659756) q[2];
sx q[2];
rz(-1.6483726) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4213688) q[1];
sx q[1];
rz(-1.6744057) q[1];
sx q[1];
rz(3.1006579) q[1];
rz(-0.85785474) q[3];
sx q[3];
rz(-0.80227533) q[3];
sx q[3];
rz(-1.0209393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3516922) q[2];
sx q[2];
rz(-1.60195) q[2];
sx q[2];
rz(-0.44463012) q[2];
rz(0.45186684) q[3];
sx q[3];
rz(-0.63932747) q[3];
sx q[3];
rz(-3.0795081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43279466) q[0];
sx q[0];
rz(-0.72467703) q[0];
sx q[0];
rz(-0.92025796) q[0];
rz(0.96011773) q[1];
sx q[1];
rz(-1.7476387) q[1];
sx q[1];
rz(0.49930176) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47682522) q[0];
sx q[0];
rz(-1.9070325) q[0];
sx q[0];
rz(1.7306843) q[0];
rz(-pi) q[1];
rz(0.22412207) q[2];
sx q[2];
rz(-2.4680228) q[2];
sx q[2];
rz(-2.3526255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.074452049) q[1];
sx q[1];
rz(-2.4765091) q[1];
sx q[1];
rz(-2.6159899) q[1];
x q[2];
rz(-1.9176513) q[3];
sx q[3];
rz(-2.2815274) q[3];
sx q[3];
rz(-0.06547346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7509191) q[2];
sx q[2];
rz(-1.5385188) q[2];
sx q[2];
rz(1.0323367) q[2];
rz(2.27683) q[3];
sx q[3];
rz(-1.4356177) q[3];
sx q[3];
rz(2.5800956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.7086577) q[0];
sx q[0];
rz(-1.9517887) q[0];
sx q[0];
rz(-1.7565961) q[0];
rz(2.2191091) q[1];
sx q[1];
rz(-1.9360767) q[1];
sx q[1];
rz(-0.14499697) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50853106) q[0];
sx q[0];
rz(-0.2446188) q[0];
sx q[0];
rz(-0.82392719) q[0];
rz(-0.85741557) q[2];
sx q[2];
rz(-1.5930297) q[2];
sx q[2];
rz(-2.2043101) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0516982) q[1];
sx q[1];
rz(-1.5400008) q[1];
sx q[1];
rz(-1.4975862) q[1];
rz(-pi) q[2];
rz(2.1337521) q[3];
sx q[3];
rz(-2.2641695) q[3];
sx q[3];
rz(-1.0348606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21236913) q[2];
sx q[2];
rz(-2.9085458) q[2];
sx q[2];
rz(1.2082427) q[2];
rz(-0.82331795) q[3];
sx q[3];
rz(-1.1813141) q[3];
sx q[3];
rz(-1.7821144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.061148297) q[0];
sx q[0];
rz(-0.7875945) q[0];
sx q[0];
rz(1.0795235) q[0];
rz(-2.1615255) q[1];
sx q[1];
rz(-1.020224) q[1];
sx q[1];
rz(-3.0316839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85724505) q[0];
sx q[0];
rz(-1.9316512) q[0];
sx q[0];
rz(-1.2566823) q[0];
x q[1];
rz(0.27836043) q[2];
sx q[2];
rz(-1.5794069) q[2];
sx q[2];
rz(1.4573163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3617864) q[1];
sx q[1];
rz(-2.61269) q[1];
sx q[1];
rz(2.5216549) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1004543) q[3];
sx q[3];
rz(-1.7653437) q[3];
sx q[3];
rz(1.1916849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2073233) q[2];
sx q[2];
rz(-1.0730275) q[2];
sx q[2];
rz(-0.70460021) q[2];
rz(-2.8568824) q[3];
sx q[3];
rz(-2.4094818) q[3];
sx q[3];
rz(2.708191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1339742) q[0];
sx q[0];
rz(-1.824279) q[0];
sx q[0];
rz(0.91868573) q[0];
rz(2.6559415) q[1];
sx q[1];
rz(-0.44356569) q[1];
sx q[1];
rz(-2.0408911) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2393409) q[0];
sx q[0];
rz(-1.709769) q[0];
sx q[0];
rz(-3.001717) q[0];
rz(0.23156872) q[2];
sx q[2];
rz(-1.8652038) q[2];
sx q[2];
rz(-0.60588479) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6599258) q[1];
sx q[1];
rz(-1.9710014) q[1];
sx q[1];
rz(3.0277962) q[1];
rz(-pi) q[2];
rz(-1.6343377) q[3];
sx q[3];
rz(-0.93616928) q[3];
sx q[3];
rz(-1.612048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6653768) q[2];
sx q[2];
rz(-2.2874338) q[2];
sx q[2];
rz(0.7594792) q[2];
rz(2.1935943) q[3];
sx q[3];
rz(-1.6839323) q[3];
sx q[3];
rz(-1.6215526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043902472) q[0];
sx q[0];
rz(-0.49711415) q[0];
sx q[0];
rz(-0.33101606) q[0];
rz(0.76200062) q[1];
sx q[1];
rz(-2.8490729) q[1];
sx q[1];
rz(-0.62823546) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10234235) q[0];
sx q[0];
rz(-1.2919173) q[0];
sx q[0];
rz(2.5526702) q[0];
rz(-pi) q[1];
rz(-2.4961595) q[2];
sx q[2];
rz(-2.7629768) q[2];
sx q[2];
rz(-2.9032674) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8208295) q[1];
sx q[1];
rz(-0.71457982) q[1];
sx q[1];
rz(0.54820593) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9167494) q[3];
sx q[3];
rz(-1.9613492) q[3];
sx q[3];
rz(-1.8486763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1198173) q[2];
sx q[2];
rz(-0.47406083) q[2];
sx q[2];
rz(-0.21381703) q[2];
rz(-0.98624936) q[3];
sx q[3];
rz(-1.8503559) q[3];
sx q[3];
rz(-0.10828644) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66961359) q[0];
sx q[0];
rz(-1.6983953) q[0];
sx q[0];
rz(-1.4629913) q[0];
rz(1.1769453) q[1];
sx q[1];
rz(-1.6564449) q[1];
sx q[1];
rz(0.0080531837) q[1];
rz(0.17433737) q[2];
sx q[2];
rz(-1.2703036) q[2];
sx q[2];
rz(-1.6008137) q[2];
rz(-1.6515274) q[3];
sx q[3];
rz(-2.2207618) q[3];
sx q[3];
rz(-0.13151463) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
