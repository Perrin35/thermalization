OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(1.6605641) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(3.0552157) q[1];
sx q[1];
rz(9.4019158) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4532783) q[0];
sx q[0];
rz(-0.18916057) q[0];
sx q[0];
rz(2.2384089) q[0];
rz(0.68732287) q[2];
sx q[2];
rz(-2.5097373) q[2];
sx q[2];
rz(-0.99422115) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37576807) q[1];
sx q[1];
rz(-1.4210912) q[1];
sx q[1];
rz(-1.9181812) q[1];
x q[2];
rz(-2.1894737) q[3];
sx q[3];
rz(-1.3631571) q[3];
sx q[3];
rz(-1.4116956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(0.76618761) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3294753) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(-2.2791729) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(-0.08509732) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9239685) q[0];
sx q[0];
rz(-2.1930165) q[0];
sx q[0];
rz(-1.5505962) q[0];
rz(-pi) q[1];
rz(0.30479635) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(-3.0457029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9566006) q[1];
sx q[1];
rz(-1.6829832) q[1];
sx q[1];
rz(0.086602028) q[1];
rz(-pi) q[2];
rz(2.8220903) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(0.055469661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8791085) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(1.6764199) q[2];
rz(1.9654467) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(-2.235967) q[0];
rz(0.33272818) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.3844301) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8850088) q[0];
sx q[0];
rz(-0.48261595) q[0];
sx q[0];
rz(-2.9628739) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7483814) q[2];
sx q[2];
rz(-2.0438497) q[2];
sx q[2];
rz(-2.2676603) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78450173) q[1];
sx q[1];
rz(-0.75575268) q[1];
sx q[1];
rz(-0.032553629) q[1];
rz(-pi) q[2];
rz(-1.0043886) q[3];
sx q[3];
rz(-1.1096138) q[3];
sx q[3];
rz(-3.0571292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(-1.8133694) q[2];
rz(1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(-3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06552799) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(2.4705825) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(0.20733325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3509388) q[0];
sx q[0];
rz(-1.6121943) q[0];
sx q[0];
rz(-3.1326276) q[0];
x q[1];
rz(-1.3968799) q[2];
sx q[2];
rz(-1.6007649) q[2];
sx q[2];
rz(-2.7001691) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3934717) q[1];
sx q[1];
rz(-1.3337743) q[1];
sx q[1];
rz(-2.9841828) q[1];
rz(-pi) q[2];
rz(1.5563577) q[3];
sx q[3];
rz(-2.0824021) q[3];
sx q[3];
rz(2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8308782) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(-0.26516178) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(2.6079544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424608) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(1.2752227) q[0];
rz(-pi) q[1];
rz(-2.3737714) q[2];
sx q[2];
rz(-1.4066833) q[2];
sx q[2];
rz(-2.5196911) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4662019) q[1];
sx q[1];
rz(-1.7168105) q[1];
sx q[1];
rz(-2.9195665) q[1];
rz(2.93612) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(-1.458414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-0.56838244) q[2];
rz(-1.9249632) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(-1.0046545) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6509103) q[0];
sx q[0];
rz(-2.6255529) q[0];
sx q[0];
rz(0.74923058) q[0];
rz(-0.3968233) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(-0.65149388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.172563) q[1];
sx q[1];
rz(-1.3247715) q[1];
sx q[1];
rz(0.56117705) q[1];
rz(-pi) q[2];
rz(1.1569571) q[3];
sx q[3];
rz(-1.5817747) q[3];
sx q[3];
rz(-1.4423086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7375609) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77666831) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(0.65814322) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(2.4694494) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015398) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(0.0075017651) q[0];
rz(-pi) q[1];
rz(-0.80588801) q[2];
sx q[2];
rz(-0.58008367) q[2];
sx q[2];
rz(-2.9758331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0828515) q[1];
sx q[1];
rz(-2.5141659) q[1];
sx q[1];
rz(0.65545603) q[1];
rz(1.4268731) q[3];
sx q[3];
rz(-0.87464273) q[3];
sx q[3];
rz(1.6203984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5333574) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(-2.3170025) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(-0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010715) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(-1.0868866) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(2.3349082) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54803941) q[0];
sx q[0];
rz(-1.1386477) q[0];
sx q[0];
rz(-3.0120864) q[0];
rz(-1.6786472) q[2];
sx q[2];
rz(-1.7753522) q[2];
sx q[2];
rz(-2.4754935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1572691) q[1];
sx q[1];
rz(-1.2829797) q[1];
sx q[1];
rz(0.70043522) q[1];
x q[2];
rz(0.032012352) q[3];
sx q[3];
rz(-1.8714928) q[3];
sx q[3];
rz(2.7923982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(-0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(0.34067571) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1239615) q[0];
sx q[0];
rz(-0.88273662) q[0];
sx q[0];
rz(-0.4171564) q[0];
x q[1];
rz(1.0057783) q[2];
sx q[2];
rz(-1.7592332) q[2];
sx q[2];
rz(0.68067683) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9790736) q[1];
sx q[1];
rz(-2.5538429) q[1];
sx q[1];
rz(-1.2390922) q[1];
rz(-pi) q[2];
rz(-2.4386028) q[3];
sx q[3];
rz(-0.25745108) q[3];
sx q[3];
rz(-1.6182181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(-0.49794751) q[2];
rz(2.8399816) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-0.27858946) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3007606) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(-1.8873909) q[0];
x q[1];
rz(-2.9701091) q[2];
sx q[2];
rz(-1.0043317) q[2];
sx q[2];
rz(-0.4085853) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6381166) q[1];
sx q[1];
rz(-1.5498811) q[1];
sx q[1];
rz(0.36550826) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46080188) q[3];
sx q[3];
rz(-2.9446903) q[3];
sx q[3];
rz(-0.27289665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(2.0417487) q[2];
sx q[2];
rz(-1.8269314) q[2];
sx q[2];
rz(-1.6614428) q[2];
rz(1.9361817) q[3];
sx q[3];
rz(-1.1689417) q[3];
sx q[3];
rz(0.4369215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
