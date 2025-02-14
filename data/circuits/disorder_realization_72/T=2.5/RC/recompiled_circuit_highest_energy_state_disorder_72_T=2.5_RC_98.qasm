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
rz(0.82062757) q[0];
sx q[0];
rz(2.4790915) q[0];
sx q[0];
rz(9.6883246) q[0];
rz(1.1072371) q[1];
sx q[1];
rz(-1.8594445) q[1];
sx q[1];
rz(0.67556226) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0209262) q[0];
sx q[0];
rz(-1.5643189) q[0];
sx q[0];
rz(-1.2351888) q[0];
rz(-pi) q[1];
rz(-2.1525429) q[2];
sx q[2];
rz(-2.0768696) q[2];
sx q[2];
rz(-0.29935867) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0274715) q[1];
sx q[1];
rz(-2.7222789) q[1];
sx q[1];
rz(-1.9472576) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3278424) q[3];
sx q[3];
rz(-1.0910506) q[3];
sx q[3];
rz(-2.858851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2519309) q[2];
sx q[2];
rz(-1.9193005) q[2];
sx q[2];
rz(-0.3678073) q[2];
rz(-2.8273072) q[3];
sx q[3];
rz(-2.8511484) q[3];
sx q[3];
rz(-0.17816003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3322068) q[0];
sx q[0];
rz(-3.0520913) q[0];
sx q[0];
rz(-1.3595164) q[0];
rz(1.8571732) q[1];
sx q[1];
rz(-1.9764683) q[1];
sx q[1];
rz(3.0899835) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29362677) q[0];
sx q[0];
rz(-0.32362263) q[0];
sx q[0];
rz(-2.1725562) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.935595) q[2];
sx q[2];
rz(-1.2075181) q[2];
sx q[2];
rz(1.3000803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7588531) q[1];
sx q[1];
rz(-1.570228) q[1];
sx q[1];
rz(-1.9447536) q[1];
rz(-pi) q[2];
rz(0.47447954) q[3];
sx q[3];
rz(-1.9360376) q[3];
sx q[3];
rz(0.38680916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0525558) q[2];
sx q[2];
rz(-2.8812228) q[2];
sx q[2];
rz(0.54711771) q[2];
rz(-1.8735006) q[3];
sx q[3];
rz(-1.3601466) q[3];
sx q[3];
rz(0.33856302) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7608305) q[0];
sx q[0];
rz(-1.0365423) q[0];
sx q[0];
rz(-1.8007675) q[0];
rz(-2.2876168) q[1];
sx q[1];
rz(-1.8546591) q[1];
sx q[1];
rz(0.30240789) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15203619) q[0];
sx q[0];
rz(-1.2947417) q[0];
sx q[0];
rz(0.92854519) q[0];
x q[1];
rz(-1.339389) q[2];
sx q[2];
rz(-2.0112027) q[2];
sx q[2];
rz(-0.4612067) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93487996) q[1];
sx q[1];
rz(-1.0374962) q[1];
sx q[1];
rz(-0.7278022) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0875999) q[3];
sx q[3];
rz(-2.5977166) q[3];
sx q[3];
rz(-0.15009201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2053908) q[2];
sx q[2];
rz(-1.8855636) q[2];
sx q[2];
rz(-2.5993247) q[2];
rz(-2.1451456) q[3];
sx q[3];
rz(-2.9139329) q[3];
sx q[3];
rz(0.12797932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.19313136) q[0];
sx q[0];
rz(-1.6574991) q[0];
sx q[0];
rz(1.2633854) q[0];
rz(-2.7073233) q[1];
sx q[1];
rz(-0.77249211) q[1];
sx q[1];
rz(-1.6901406) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5853441) q[0];
sx q[0];
rz(-0.55216575) q[0];
sx q[0];
rz(0.57741965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66777467) q[2];
sx q[2];
rz(-0.66168565) q[2];
sx q[2];
rz(-1.1061321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60148224) q[1];
sx q[1];
rz(-0.14232351) q[1];
sx q[1];
rz(2.2156634) q[1];
rz(-pi) q[2];
rz(-0.94971117) q[3];
sx q[3];
rz(-1.1009645) q[3];
sx q[3];
rz(-2.7787152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0749391) q[2];
sx q[2];
rz(-2.1982919) q[2];
sx q[2];
rz(2.2554876) q[2];
rz(0.0063627176) q[3];
sx q[3];
rz(-1.3402091) q[3];
sx q[3];
rz(1.1191012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5804382) q[0];
sx q[0];
rz(-2.6851324) q[0];
sx q[0];
rz(-1.891267) q[0];
rz(2.9087032) q[1];
sx q[1];
rz(-0.75585514) q[1];
sx q[1];
rz(-3.0511391) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2523777) q[0];
sx q[0];
rz(-1.4486635) q[0];
sx q[0];
rz(-1.8218763) q[0];
rz(2.6695741) q[2];
sx q[2];
rz(-1.8119805) q[2];
sx q[2];
rz(-2.3015353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.662218) q[1];
sx q[1];
rz(-1.5081128) q[1];
sx q[1];
rz(-2.4880243) q[1];
x q[2];
rz(-1.8985604) q[3];
sx q[3];
rz(-1.6521427) q[3];
sx q[3];
rz(1.3535172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6470486) q[2];
sx q[2];
rz(-2.9477951) q[2];
sx q[2];
rz(-1.2079283) q[2];
rz(-2.2203994) q[3];
sx q[3];
rz(-2.2974206) q[3];
sx q[3];
rz(0.45879656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0334466) q[0];
sx q[0];
rz(-1.7195846) q[0];
sx q[0];
rz(2.6500927) q[0];
rz(0.72832251) q[1];
sx q[1];
rz(-1.1110543) q[1];
sx q[1];
rz(-0.19457766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.174343) q[0];
sx q[0];
rz(-1.3493269) q[0];
sx q[0];
rz(-2.1465149) q[0];
rz(2.3512882) q[2];
sx q[2];
rz(-0.93043295) q[2];
sx q[2];
rz(1.9483639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.814333) q[1];
sx q[1];
rz(-0.93686283) q[1];
sx q[1];
rz(2.4933706) q[1];
rz(-1.3380906) q[3];
sx q[3];
rz(-0.92706087) q[3];
sx q[3];
rz(0.88730592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7877833) q[2];
sx q[2];
rz(-2.1259191) q[2];
sx q[2];
rz(3.1397528) q[2];
rz(1.4932102) q[3];
sx q[3];
rz(-0.73072481) q[3];
sx q[3];
rz(-0.28436896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581185) q[0];
sx q[0];
rz(-2.879877) q[0];
sx q[0];
rz(2.0765685) q[0];
rz(2.8788772) q[1];
sx q[1];
rz(-2.5698667) q[1];
sx q[1];
rz(0.46331847) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36998591) q[0];
sx q[0];
rz(-0.60877234) q[0];
sx q[0];
rz(1.8474402) q[0];
rz(-2.1049785) q[2];
sx q[2];
rz(-2.2232703) q[2];
sx q[2];
rz(-2.5135615) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66760705) q[1];
sx q[1];
rz(-1.7069478) q[1];
sx q[1];
rz(2.6546247) q[1];
x q[2];
rz(1.7780532) q[3];
sx q[3];
rz(-0.1175783) q[3];
sx q[3];
rz(0.97299796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.00039438417) q[2];
sx q[2];
rz(-1.1606777) q[2];
sx q[2];
rz(2.311643) q[2];
rz(1.0478033) q[3];
sx q[3];
rz(-0.44568291) q[3];
sx q[3];
rz(-2.9629663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038427) q[0];
sx q[0];
rz(-0.99030817) q[0];
sx q[0];
rz(2.3867699) q[0];
rz(-2.7524475) q[1];
sx q[1];
rz(-1.6837589) q[1];
sx q[1];
rz(2.7438502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75466418) q[0];
sx q[0];
rz(-0.3873741) q[0];
sx q[0];
rz(-2.3088916) q[0];
rz(-1.7988458) q[2];
sx q[2];
rz(-1.031165) q[2];
sx q[2];
rz(-1.9613105) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6036247) q[1];
sx q[1];
rz(-1.84314) q[1];
sx q[1];
rz(-2.8773099) q[1];
rz(-pi) q[2];
rz(-0.92614545) q[3];
sx q[3];
rz(-2.20813) q[3];
sx q[3];
rz(1.6660575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3079181) q[2];
sx q[2];
rz(-2.8454744) q[2];
sx q[2];
rz(-2.8274242) q[2];
rz(1.2416154) q[3];
sx q[3];
rz(-1.5528468) q[3];
sx q[3];
rz(-1.9023021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32387963) q[0];
sx q[0];
rz(-2.4535024) q[0];
sx q[0];
rz(-2.896198) q[0];
rz(-1.2303526) q[1];
sx q[1];
rz(-0.10086682) q[1];
sx q[1];
rz(2.6255677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8300872) q[0];
sx q[0];
rz(-1.9507512) q[0];
sx q[0];
rz(1.3561234) q[0];
rz(-pi) q[1];
rz(-2.3886959) q[2];
sx q[2];
rz(-1.0123526) q[2];
sx q[2];
rz(3.0679997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5616331) q[1];
sx q[1];
rz(-1.206534) q[1];
sx q[1];
rz(1.2666525) q[1];
x q[2];
rz(3.0050659) q[3];
sx q[3];
rz(-1.6322005) q[3];
sx q[3];
rz(-1.2207954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5838991) q[2];
sx q[2];
rz(-1.9198753) q[2];
sx q[2];
rz(-2.058513) q[2];
rz(0.22100581) q[3];
sx q[3];
rz(-0.14328863) q[3];
sx q[3];
rz(-0.7964645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021127064) q[0];
sx q[0];
rz(-3.0607304) q[0];
sx q[0];
rz(3.0057111) q[0];
rz(1.7120301) q[1];
sx q[1];
rz(-2.0202961) q[1];
sx q[1];
rz(2.8543465) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82086876) q[0];
sx q[0];
rz(-2.831376) q[0];
sx q[0];
rz(2.2287772) q[0];
rz(3.0570266) q[2];
sx q[2];
rz(-1.8075602) q[2];
sx q[2];
rz(2.6344476) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96480709) q[1];
sx q[1];
rz(-0.63399708) q[1];
sx q[1];
rz(0.37668677) q[1];
rz(-pi) q[2];
rz(0.29218896) q[3];
sx q[3];
rz(-1.0064408) q[3];
sx q[3];
rz(2.5404504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2222432) q[2];
sx q[2];
rz(-1.9025849) q[2];
sx q[2];
rz(-2.3075721) q[2];
rz(-0.30759865) q[3];
sx q[3];
rz(-2.7643876) q[3];
sx q[3];
rz(-0.26267499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.96497) q[0];
sx q[0];
rz(-1.3055834) q[0];
sx q[0];
rz(2.1156043) q[0];
rz(2.7550244) q[1];
sx q[1];
rz(-1.3810806) q[1];
sx q[1];
rz(2.5234533) q[1];
rz(1.7025399) q[2];
sx q[2];
rz(-1.1494888) q[2];
sx q[2];
rz(-0.26486808) q[2];
rz(-1.3597506) q[3];
sx q[3];
rz(-1.4283709) q[3];
sx q[3];
rz(-1.4621468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
