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
rz(-2.8060198) q[0];
sx q[0];
rz(-1.6874474) q[0];
sx q[0];
rz(-1.0336706) q[0];
rz(-1.7103065) q[1];
sx q[1];
rz(-2.5841076) q[1];
sx q[1];
rz(-2.1318336) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.570734) q[0];
sx q[0];
rz(-1.6964127) q[0];
sx q[0];
rz(-1.3799589) q[0];
rz(-pi) q[1];
rz(1.5255341) q[2];
sx q[2];
rz(-0.082556225) q[2];
sx q[2];
rz(-2.1545067) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.941135) q[1];
sx q[1];
rz(-1.7361281) q[1];
sx q[1];
rz(-1.6881866) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5663773) q[3];
sx q[3];
rz(-1.7130245) q[3];
sx q[3];
rz(2.9517236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7734163) q[2];
sx q[2];
rz(-0.877031) q[2];
sx q[2];
rz(-2.3167493) q[2];
rz(-0.33251897) q[3];
sx q[3];
rz(-0.85086346) q[3];
sx q[3];
rz(-0.48377812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2291718) q[0];
sx q[0];
rz(-0.084675463) q[0];
sx q[0];
rz(0.53572768) q[0];
rz(2.3748705) q[1];
sx q[1];
rz(-1.352939) q[1];
sx q[1];
rz(1.3923233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6321997) q[0];
sx q[0];
rz(-2.0344816) q[0];
sx q[0];
rz(-0.38393387) q[0];
x q[1];
rz(-3.1145623) q[2];
sx q[2];
rz(-0.34124869) q[2];
sx q[2];
rz(0.42021593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.036511985) q[1];
sx q[1];
rz(-2.3859508) q[1];
sx q[1];
rz(-1.6022268) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6757319) q[3];
sx q[3];
rz(-1.6989048) q[3];
sx q[3];
rz(2.6554633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2367737) q[2];
sx q[2];
rz(-2.3084013) q[2];
sx q[2];
rz(-2.0429677) q[2];
rz(2.3084579) q[3];
sx q[3];
rz(-1.5980709) q[3];
sx q[3];
rz(-1.0540849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5611834) q[0];
sx q[0];
rz(-1.9920749) q[0];
sx q[0];
rz(2.4533601) q[0];
rz(-1.8904103) q[1];
sx q[1];
rz(-2.4871608) q[1];
sx q[1];
rz(-1.2122663) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040846856) q[0];
sx q[0];
rz(-1.0991352) q[0];
sx q[0];
rz(-0.77420401) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86264865) q[2];
sx q[2];
rz(-2.7852163) q[2];
sx q[2];
rz(-1.4087806) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40078354) q[1];
sx q[1];
rz(-0.91183981) q[1];
sx q[1];
rz(-1.5636921) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5098677) q[3];
sx q[3];
rz(-3.0620771) q[3];
sx q[3];
rz(2.7729976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4091829) q[2];
sx q[2];
rz(-2.1467291) q[2];
sx q[2];
rz(-0.44962064) q[2];
rz(-2.2879587) q[3];
sx q[3];
rz(-0.64695224) q[3];
sx q[3];
rz(0.70931119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.7483826) q[0];
sx q[0];
rz(-2.125183) q[0];
sx q[0];
rz(2.7384695) q[0];
rz(-2.368811) q[1];
sx q[1];
rz(-1.9085596) q[1];
sx q[1];
rz(-2.3162139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7220424) q[0];
sx q[0];
rz(-0.78762509) q[0];
sx q[0];
rz(-1.1088637) q[0];
rz(-3.1146997) q[2];
sx q[2];
rz(-1.5468569) q[2];
sx q[2];
rz(2.8223512) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1267438) q[1];
sx q[1];
rz(-0.9909317) q[1];
sx q[1];
rz(-0.090437263) q[1];
rz(-2.585936) q[3];
sx q[3];
rz(-0.98816493) q[3];
sx q[3];
rz(2.5224796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0394502) q[2];
sx q[2];
rz(-1.2692229) q[2];
sx q[2];
rz(-0.5086745) q[2];
rz(-1.8968808) q[3];
sx q[3];
rz(-2.0679943) q[3];
sx q[3];
rz(-1.8083474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4077069) q[0];
sx q[0];
rz(-1.5835967) q[0];
sx q[0];
rz(-0.35633126) q[0];
rz(-1.7286667) q[1];
sx q[1];
rz(-1.3099542) q[1];
sx q[1];
rz(-1.2948571) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0771478) q[0];
sx q[0];
rz(-1.1379717) q[0];
sx q[0];
rz(-1.8124142) q[0];
x q[1];
rz(0.21759502) q[2];
sx q[2];
rz(-2.4684973) q[2];
sx q[2];
rz(-0.714528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8581681) q[1];
sx q[1];
rz(-1.8008261) q[1];
sx q[1];
rz(-2.7344804) q[1];
rz(-pi) q[2];
rz(-1.0255085) q[3];
sx q[3];
rz(-1.9615615) q[3];
sx q[3];
rz(2.7440967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9406708) q[2];
sx q[2];
rz(-1.385043) q[2];
sx q[2];
rz(0.6400288) q[2];
rz(-2.707543) q[3];
sx q[3];
rz(-0.95165747) q[3];
sx q[3];
rz(-1.1205589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5956748) q[0];
sx q[0];
rz(-2.7037103) q[0];
sx q[0];
rz(-0.75089279) q[0];
rz(-0.76260507) q[1];
sx q[1];
rz(-1.7060988) q[1];
sx q[1];
rz(2.6079752) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.740858) q[0];
sx q[0];
rz(-2.7228231) q[0];
sx q[0];
rz(-0.77218582) q[0];
rz(0.80383382) q[2];
sx q[2];
rz(-1.8920004) q[2];
sx q[2];
rz(0.28576947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4009559) q[1];
sx q[1];
rz(-1.4527241) q[1];
sx q[1];
rz(-0.90098874) q[1];
rz(2.153965) q[3];
sx q[3];
rz(-1.0915874) q[3];
sx q[3];
rz(-0.85239172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.077347191) q[2];
sx q[2];
rz(-2.7392445) q[2];
sx q[2];
rz(0.8858436) q[2];
rz(-1.8387851) q[3];
sx q[3];
rz(-1.9582483) q[3];
sx q[3];
rz(2.6634789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2122022) q[0];
sx q[0];
rz(-2.2181856) q[0];
sx q[0];
rz(2.5481664) q[0];
rz(1.5133096) q[1];
sx q[1];
rz(-1.9624886) q[1];
sx q[1];
rz(2.5679307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5710756) q[0];
sx q[0];
rz(-0.62643753) q[0];
sx q[0];
rz(2.4721739) q[0];
rz(-pi) q[1];
rz(-1.8183579) q[2];
sx q[2];
rz(-1.0070966) q[2];
sx q[2];
rz(-0.42681387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5630514) q[1];
sx q[1];
rz(-2.6669901) q[1];
sx q[1];
rz(-2.5792349) q[1];
rz(1.6214293) q[3];
sx q[3];
rz(-1.4922499) q[3];
sx q[3];
rz(-2.1861064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71184984) q[2];
sx q[2];
rz(-1.9387127) q[2];
sx q[2];
rz(1.4329866) q[2];
rz(1.1687219) q[3];
sx q[3];
rz(-2.70372) q[3];
sx q[3];
rz(-1.6247113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.8716005) q[0];
sx q[0];
rz(-0.90684909) q[0];
sx q[0];
rz(-2.9378743) q[0];
rz(-1.3153971) q[1];
sx q[1];
rz(-1.8351277) q[1];
sx q[1];
rz(-1.809583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4341605) q[0];
sx q[0];
rz(-0.95904175) q[0];
sx q[0];
rz(0.79025288) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88896759) q[2];
sx q[2];
rz(-2.3485514) q[2];
sx q[2];
rz(-0.5785795) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9086919) q[1];
sx q[1];
rz(-0.76592839) q[1];
sx q[1];
rz(0.58118622) q[1];
x q[2];
rz(-0.90059728) q[3];
sx q[3];
rz(-0.20684563) q[3];
sx q[3];
rz(-1.2634009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3709162) q[2];
sx q[2];
rz(-0.081826536) q[2];
sx q[2];
rz(-1.3814242) q[2];
rz(-0.55111876) q[3];
sx q[3];
rz(-1.8128606) q[3];
sx q[3];
rz(2.3962925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4296752) q[0];
sx q[0];
rz(-1.4497919) q[0];
sx q[0];
rz(1.2592738) q[0];
rz(-1.7970386) q[1];
sx q[1];
rz(-0.63251907) q[1];
sx q[1];
rz(-1.22619) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38905242) q[0];
sx q[0];
rz(-1.2864094) q[0];
sx q[0];
rz(0.14099462) q[0];
x q[1];
rz(-1.1414503) q[2];
sx q[2];
rz(-1.2163905) q[2];
sx q[2];
rz(2.5130075) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.256168) q[1];
sx q[1];
rz(-2.3162335) q[1];
sx q[1];
rz(-1.481545) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1192083) q[3];
sx q[3];
rz(-2.5520241) q[3];
sx q[3];
rz(-2.0062674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5598477) q[2];
sx q[2];
rz(-0.86220828) q[2];
sx q[2];
rz(-1.6017412) q[2];
rz(-1.3567443) q[3];
sx q[3];
rz(-0.53250766) q[3];
sx q[3];
rz(1.2342359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27720472) q[0];
sx q[0];
rz(-1.5437523) q[0];
sx q[0];
rz(-3.0371015) q[0];
rz(-1.744572) q[1];
sx q[1];
rz(-1.3518159) q[1];
sx q[1];
rz(-1.6207961) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0919017) q[0];
sx q[0];
rz(-0.77020634) q[0];
sx q[0];
rz(0.30253221) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42664032) q[2];
sx q[2];
rz(-0.37084118) q[2];
sx q[2];
rz(-1.1877354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0238259) q[1];
sx q[1];
rz(-1.858288) q[1];
sx q[1];
rz(1.4625129) q[1];
rz(-1.4827221) q[3];
sx q[3];
rz(-0.58335005) q[3];
sx q[3];
rz(-0.47881918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4349159) q[2];
sx q[2];
rz(-0.84182635) q[2];
sx q[2];
rz(-1.288877) q[2];
rz(-2.1364818) q[3];
sx q[3];
rz(-1.0594599) q[3];
sx q[3];
rz(3.03481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91916753) q[0];
sx q[0];
rz(-1.5594302) q[0];
sx q[0];
rz(2.9710309) q[0];
rz(-0.87986058) q[1];
sx q[1];
rz(-0.51645551) q[1];
sx q[1];
rz(-2.5194306) q[1];
rz(1.6352898) q[2];
sx q[2];
rz(-2.0694824) q[2];
sx q[2];
rz(2.4691442) q[2];
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
