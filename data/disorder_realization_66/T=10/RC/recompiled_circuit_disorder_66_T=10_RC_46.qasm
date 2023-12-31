OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(4.1132676) q[0];
sx q[0];
rz(10.905807) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(0.02286214) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6883144) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(-2.2384089) q[0];
rz(-pi) q[1];
rz(-2.0055662) q[2];
sx q[2];
rz(-1.0966986) q[2];
sx q[2];
rz(-2.9413162) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80427985) q[1];
sx q[1];
rz(-0.37706456) q[1];
sx q[1];
rz(-1.9878597) q[1];
x q[2];
rz(-1.919235) q[3];
sx q[3];
rz(-0.64823965) q[3];
sx q[3];
rz(3.0188308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8121174) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-0.086659327) q[0];
rz(2.2791729) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(3.0564953) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25227308) q[0];
sx q[0];
rz(-2.5190881) q[0];
sx q[0];
rz(-3.1134393) q[0];
rz(-pi) q[1];
rz(-0.30479635) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(3.0457029) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18499204) q[1];
sx q[1];
rz(-1.4586095) q[1];
sx q[1];
rz(-3.0549906) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1477633) q[3];
sx q[3];
rz(-0.54518632) q[3];
sx q[3];
rz(0.59515566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(-2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.3844301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8850088) q[0];
sx q[0];
rz(-0.48261595) q[0];
sx q[0];
rz(-2.9628739) q[0];
rz(-pi) q[1];
rz(-0.92817523) q[2];
sx q[2];
rz(-2.536142) q[2];
sx q[2];
rz(1.5290934) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82921925) q[1];
sx q[1];
rz(-2.3260498) q[1];
sx q[1];
rz(1.5401328) q[1];
rz(-pi) q[2];
rz(2.3178187) q[3];
sx q[3];
rz(-2.4274821) q[3];
sx q[3];
rz(0.87574524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.3282233) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06552799) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(-1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(2.9342594) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78051356) q[0];
sx q[0];
rz(-1.5797537) q[0];
sx q[0];
rz(1.6121959) q[0];
rz(3.1111654) q[2];
sx q[2];
rz(-1.7446339) q[2];
sx q[2];
rz(-1.1241084) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15370788) q[1];
sx q[1];
rz(-0.28370902) q[1];
sx q[1];
rz(0.99516408) q[1];
rz(-pi) q[2];
rz(-0.51165032) q[3];
sx q[3];
rz(-1.5582065) q[3];
sx q[3];
rz(0.97012855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31071445) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-0.19367735) q[3];
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
rz(0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(-0.26516178) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(2.6079544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9981209) q[0];
sx q[0];
rz(-0.40092418) q[0];
sx q[0];
rz(-0.801416) q[0];
rz(-pi) q[1];
rz(0.76782121) q[2];
sx q[2];
rz(-1.7349093) q[2];
sx q[2];
rz(2.5196911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4662019) q[1];
sx q[1];
rz(-1.4247822) q[1];
sx q[1];
rz(0.22202613) q[1];
rz(-pi) q[2];
rz(-2.93612) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(1.458414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9020033) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-0.56838244) q[2];
rz(1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(2.7980428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6333273) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(-1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(1.0046545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4906824) q[0];
sx q[0];
rz(-2.6255529) q[0];
sx q[0];
rz(-0.74923058) q[0];
rz(-pi) q[1];
rz(2.7447694) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(2.4900988) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.172563) q[1];
sx q[1];
rz(-1.3247715) q[1];
sx q[1];
rz(-0.56117705) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5435013) q[3];
sx q[3];
rz(-2.7276162) q[3];
sx q[3];
rz(-2.9881145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7375609) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(2.693434) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-2.7929849) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(-0.67214322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2015398) q[0];
sx q[0];
rz(-1.4048028) q[0];
sx q[0];
rz(-3.1340909) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7156098) q[2];
sx q[2];
rz(-1.1642712) q[2];
sx q[2];
rz(2.4533518) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8183648) q[1];
sx q[1];
rz(-2.0548901) q[1];
sx q[1];
rz(1.9869884) q[1];
rz(-pi) q[2];
rz(0.70127212) q[3];
sx q[3];
rz(-1.6810732) q[3];
sx q[3];
rz(2.9993204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5333574) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(-2.7427924) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44052112) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-0.80668443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0643443) q[0];
sx q[0];
rz(-1.6883388) q[0];
sx q[0];
rz(-2.006152) q[0];
x q[1];
rz(0.47862349) q[2];
sx q[2];
rz(-0.23089409) q[2];
sx q[2];
rz(2.965197) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.17864171) q[1];
sx q[1];
rz(-0.9045524) q[1];
sx q[1];
rz(-1.9402177) q[1];
x q[2];
rz(3.1095803) q[3];
sx q[3];
rz(-1.8714928) q[3];
sx q[3];
rz(-2.7923982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(0.99964833) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(-1.1243533) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12061159) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(-0.34067571) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4083456) q[0];
sx q[0];
rz(-0.78662965) q[0];
sx q[0];
rz(1.1128845) q[0];
x q[1];
rz(-1.9129487) q[2];
sx q[2];
rz(-0.59235448) q[2];
sx q[2];
rz(0.60281384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.454168) q[1];
sx q[1];
rz(-1.3892281) q[1];
sx q[1];
rz(-2.1329692) q[1];
rz(-pi) q[2];
rz(-2.4386028) q[3];
sx q[3];
rz(-0.25745108) q[3];
sx q[3];
rz(-1.6182181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3866773) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(2.6436451) q[2];
rz(2.8399816) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(0.27858946) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(3.0648807) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4519276) q[0];
sx q[0];
rz(-1.2566914) q[0];
sx q[0];
rz(0.12977022) q[0];
rz(-pi) q[1];
rz(-0.99761988) q[2];
sx q[2];
rz(-1.7152889) q[2];
sx q[2];
rz(1.0695374) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6381166) q[1];
sx q[1];
rz(-1.5498811) q[1];
sx q[1];
rz(-2.7760844) q[1];
x q[2];
rz(-2.9647787) q[3];
sx q[3];
rz(-1.483695) q[3];
sx q[3];
rz(-2.2967695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-2.6860766) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(-0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615622) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(2.5554399) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(2.0417487) q[2];
sx q[2];
rz(-1.8269314) q[2];
sx q[2];
rz(-1.6614428) q[2];
rz(2.4424845) q[3];
sx q[3];
rz(-0.53634488) q[3];
sx q[3];
rz(2.8041822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
