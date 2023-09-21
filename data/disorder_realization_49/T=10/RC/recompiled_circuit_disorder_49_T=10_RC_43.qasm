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
rz(4.0806169) q[0];
sx q[0];
rz(9.4299849) q[0];
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
rz(2.0931041) q[0];
sx q[0];
rz(-2.5079873) q[0];
sx q[0];
rz(2.5392169) q[0];
rz(-0.66733811) q[2];
sx q[2];
rz(-2.9138406) q[2];
sx q[2];
rz(1.8337133) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8813821) q[1];
sx q[1];
rz(-1.26169) q[1];
sx q[1];
rz(-0.15501546) q[1];
rz(-pi) q[2];
rz(-0.40640229) q[3];
sx q[3];
rz(-2.1592525) q[3];
sx q[3];
rz(-1.3960081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4253915) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(-2.544196) q[2];
rz(-1.3655837) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.8610154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269203) q[0];
sx q[0];
rz(-2.5898114) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(2.0479653) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(0.11322583) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118606) q[0];
sx q[0];
rz(-1.9959873) q[0];
sx q[0];
rz(3.0990764) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26248652) q[2];
sx q[2];
rz(-3.065735) q[2];
sx q[2];
rz(2.5469317) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37564056) q[1];
sx q[1];
rz(-1.1728219) q[1];
sx q[1];
rz(-0.44771938) q[1];
rz(-2.9037335) q[3];
sx q[3];
rz(-2.1616462) q[3];
sx q[3];
rz(-1.295134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.018628) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(-2.9193027) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6753321) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(-1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(-3.0139794) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022284431) q[0];
sx q[0];
rz(-0.89953178) q[0];
sx q[0];
rz(-0.22247252) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67655501) q[2];
sx q[2];
rz(-0.99202079) q[2];
sx q[2];
rz(0.81625953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0763921) q[1];
sx q[1];
rz(-0.25307357) q[1];
sx q[1];
rz(1.9051001) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2283986) q[3];
sx q[3];
rz(-0.6558334) q[3];
sx q[3];
rz(2.0369903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3601274) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(0.15790766) q[0];
rz(0.36610106) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(0.37240949) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2686553) q[0];
sx q[0];
rz(-2.4252709) q[0];
sx q[0];
rz(1.3852081) q[0];
rz(-0.33505586) q[2];
sx q[2];
rz(-2.1261566) q[2];
sx q[2];
rz(-0.43713883) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1540893) q[1];
sx q[1];
rz(-1.1945063) q[1];
sx q[1];
rz(2.980568) q[1];
rz(-2.1129205) q[3];
sx q[3];
rz(-0.91915932) q[3];
sx q[3];
rz(1.8478912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(-2.8692029) q[2];
rz(2.2327936) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(-2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(-2.0102665) q[0];
rz(2.5436026) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(-0.67684832) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15166053) q[0];
sx q[0];
rz(-2.6940072) q[0];
sx q[0];
rz(-0.18330343) q[0];
x q[1];
rz(2.9025181) q[2];
sx q[2];
rz(-0.74547807) q[2];
sx q[2];
rz(0.67666942) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21351335) q[1];
sx q[1];
rz(-0.25611862) q[1];
sx q[1];
rz(-3.03979) q[1];
rz(-pi) q[2];
rz(1.5951717) q[3];
sx q[3];
rz(-1.1773603) q[3];
sx q[3];
rz(0.76373053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1375492) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(-1.9963025) q[2];
rz(2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(-1.0173652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(1.6824678) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(2.2084592) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4695078) q[0];
sx q[0];
rz(-2.5818995) q[0];
sx q[0];
rz(2.1493388) q[0];
rz(2.6283693) q[2];
sx q[2];
rz(-0.66868082) q[2];
sx q[2];
rz(-0.19323397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16490368) q[1];
sx q[1];
rz(-0.36786825) q[1];
sx q[1];
rz(-1.5771754) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3429787) q[3];
sx q[3];
rz(-2.4711547) q[3];
sx q[3];
rz(-0.05712856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(-2.8549426) q[2];
rz(0.24946985) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1405623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46465835) q[0];
sx q[0];
rz(-1.6081728) q[0];
sx q[0];
rz(-1.8514762) q[0];
x q[1];
rz(-0.83058968) q[2];
sx q[2];
rz(-0.48362728) q[2];
sx q[2];
rz(0.73053503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8965473) q[1];
sx q[1];
rz(-0.42173112) q[1];
sx q[1];
rz(1.1827724) q[1];
rz(-pi) q[2];
rz(-1.3592968) q[3];
sx q[3];
rz(-0.930951) q[3];
sx q[3];
rz(2.0169472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.004185685) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(-0.10350791) q[2];
rz(-0.59182709) q[3];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(1.3423963) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(-2.7517095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33078157) q[0];
sx q[0];
rz(-2.1446052) q[0];
sx q[0];
rz(-1.1382889) q[0];
rz(-0.6673442) q[2];
sx q[2];
rz(-0.61605011) q[2];
sx q[2];
rz(1.6499856) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65615678) q[1];
sx q[1];
rz(-1.1404783) q[1];
sx q[1];
rz(2.7780611) q[1];
x q[2];
rz(-0.26384683) q[3];
sx q[3];
rz(-2.2328394) q[3];
sx q[3];
rz(0.86316934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-2.4273382) q[2];
rz(2.2438625) q[3];
sx q[3];
rz(-2.0102746) q[3];
sx q[3];
rz(2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(-2.3642448) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(0.27639595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24647507) q[0];
sx q[0];
rz(-2.7013489) q[0];
sx q[0];
rz(-2.759139) q[0];
x q[1];
rz(1.7644291) q[2];
sx q[2];
rz(-1.7948705) q[2];
sx q[2];
rz(-2.5335238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.249835) q[1];
sx q[1];
rz(-1.3502305) q[1];
sx q[1];
rz(2.9794429) q[1];
x q[2];
rz(-1.3235839) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(-1.1288527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2892264) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(-0.12750553) q[2];
rz(3.1048807) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(1.9071074) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7413095) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(-2.7888443) q[0];
rz(0.57669512) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(2.1113077) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0477662) q[0];
sx q[0];
rz(-0.86137912) q[0];
sx q[0];
rz(1.312027) q[0];
x q[1];
rz(-1.2376386) q[2];
sx q[2];
rz(-2.2292456) q[2];
sx q[2];
rz(-1.2308987) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5926338) q[1];
sx q[1];
rz(-2.0083798) q[1];
sx q[1];
rz(1.2195107) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.685931) q[3];
sx q[3];
rz(-0.66484287) q[3];
sx q[3];
rz(1.957422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(2.4009005) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108903) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(-0.66979349) q[1];
sx q[1];
rz(-2.5302946) q[1];
sx q[1];
rz(-1.5763462) q[1];
rz(0.40217051) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(2.4521811) q[3];
sx q[3];
rz(-1.1082311) q[3];
sx q[3];
rz(-2.5267596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];