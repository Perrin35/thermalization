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
rz(-0.79688537) q[0];
sx q[0];
rz(-1.4112043) q[0];
sx q[0];
rz(-2.222173) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(-0.04920955) q[1];
sx q[1];
rz(-0.71370178) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.365276) q[0];
sx q[0];
rz(-2.8284024) q[0];
sx q[0];
rz(1.0745144) q[0];
rz(-pi) q[1];
rz(-3.0786985) q[2];
sx q[2];
rz(-1.6990635) q[2];
sx q[2];
rz(2.4155282) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7138243) q[1];
sx q[1];
rz(-1.57219) q[1];
sx q[1];
rz(2.6906101) q[1];
x q[2];
rz(0.034293745) q[3];
sx q[3];
rz(-1.099713) q[3];
sx q[3];
rz(-0.33609875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0539703) q[2];
sx q[2];
rz(-0.38603187) q[2];
sx q[2];
rz(0.36960441) q[2];
rz(-1.9965648) q[3];
sx q[3];
rz(-1.4729045) q[3];
sx q[3];
rz(2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0358148) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(-0.57304397) q[0];
rz(2.0447958) q[1];
sx q[1];
rz(-1.6005452) q[1];
sx q[1];
rz(-2.6699452) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9100138) q[0];
sx q[0];
rz(-2.1699804) q[0];
sx q[0];
rz(3.0082361) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71895968) q[2];
sx q[2];
rz(-0.83375726) q[2];
sx q[2];
rz(1.7938839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0309033) q[1];
sx q[1];
rz(-2.2784) q[1];
sx q[1];
rz(-2.3861107) q[1];
x q[2];
rz(-2.2120958) q[3];
sx q[3];
rz(-2.5166582) q[3];
sx q[3];
rz(2.9004824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7219217) q[2];
sx q[2];
rz(-2.2025043) q[2];
sx q[2];
rz(-1.088885) q[2];
rz(-1.0645083) q[3];
sx q[3];
rz(-2.0239315) q[3];
sx q[3];
rz(2.815912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0629145) q[0];
sx q[0];
rz(-0.18904541) q[0];
sx q[0];
rz(-3.1245226) q[0];
rz(0.71006376) q[1];
sx q[1];
rz(-2.3080669) q[1];
sx q[1];
rz(2.5753218) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0294993) q[0];
sx q[0];
rz(-1.5819966) q[0];
sx q[0];
rz(2.311934) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60795075) q[2];
sx q[2];
rz(-2.3312116) q[2];
sx q[2];
rz(-1.2255526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33070811) q[1];
sx q[1];
rz(-0.20527923) q[1];
sx q[1];
rz(-0.42930557) q[1];
rz(-pi) q[2];
rz(0.2416824) q[3];
sx q[3];
rz(-0.2727601) q[3];
sx q[3];
rz(-2.5916416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8038586) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(3.1324978) q[2];
rz(-0.13633063) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(1.2135308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71753865) q[0];
sx q[0];
rz(-2.461705) q[0];
sx q[0];
rz(-0.16615443) q[0];
rz(-2.0280139) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(-0.16214935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55570468) q[0];
sx q[0];
rz(-1.7506208) q[0];
sx q[0];
rz(-0.094688133) q[0];
rz(2.5926431) q[2];
sx q[2];
rz(-1.5497583) q[2];
sx q[2];
rz(-0.32914135) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2723586) q[1];
sx q[1];
rz(-1.9290392) q[1];
sx q[1];
rz(2.3754602) q[1];
rz(-pi) q[2];
rz(2.6147551) q[3];
sx q[3];
rz(-0.73606811) q[3];
sx q[3];
rz(-1.7685304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1589511) q[2];
sx q[2];
rz(-2.0621767) q[2];
sx q[2];
rz(2.6137433) q[2];
rz(-0.85150254) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(-0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.270179) q[0];
sx q[0];
rz(-1.542792) q[0];
sx q[0];
rz(-2.340509) q[0];
rz(0.78394765) q[1];
sx q[1];
rz(-2.3990217) q[1];
sx q[1];
rz(0.98091006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1005412) q[0];
sx q[0];
rz(-1.1642191) q[0];
sx q[0];
rz(-1.8009787) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7655914) q[2];
sx q[2];
rz(-2.0472976) q[2];
sx q[2];
rz(2.8636572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6740711) q[1];
sx q[1];
rz(-1.5115807) q[1];
sx q[1];
rz(-1.6644415) q[1];
x q[2];
rz(1.1889691) q[3];
sx q[3];
rz(-2.5562048) q[3];
sx q[3];
rz(-0.25661925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.058978) q[2];
sx q[2];
rz(-1.5951944) q[2];
sx q[2];
rz(0.8832461) q[2];
rz(0.20032459) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(1.0909572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9969295) q[0];
sx q[0];
rz(-2.0893593) q[0];
sx q[0];
rz(0.99217478) q[0];
rz(0.80157533) q[1];
sx q[1];
rz(-1.8482607) q[1];
sx q[1];
rz(1.4206402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1039889) q[0];
sx q[0];
rz(-0.9524571) q[0];
sx q[0];
rz(1.9644587) q[0];
rz(3.0367766) q[2];
sx q[2];
rz(-2.1007406) q[2];
sx q[2];
rz(-1.3111834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9538624) q[1];
sx q[1];
rz(-1.8988601) q[1];
sx q[1];
rz(0.66098722) q[1];
x q[2];
rz(-1.9377557) q[3];
sx q[3];
rz(-0.10408574) q[3];
sx q[3];
rz(0.48818016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8617323) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(2.6045065) q[2];
rz(0.01072695) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(-0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2456197) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(2.8017092) q[0];
rz(1.3019568) q[1];
sx q[1];
rz(-0.94568959) q[1];
sx q[1];
rz(1.1713015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1780258) q[0];
sx q[0];
rz(-1.2739355) q[0];
sx q[0];
rz(-1.6087449) q[0];
rz(-pi) q[1];
rz(2.3717959) q[2];
sx q[2];
rz(-0.30210051) q[2];
sx q[2];
rz(-0.090902791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4170467) q[1];
sx q[1];
rz(-1.5929211) q[1];
sx q[1];
rz(1.6032277) q[1];
rz(-0.62452353) q[3];
sx q[3];
rz(-1.7134662) q[3];
sx q[3];
rz(-0.7398015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0168212) q[2];
sx q[2];
rz(-1.7243959) q[2];
sx q[2];
rz(2.4701414) q[2];
rz(-0.5365544) q[3];
sx q[3];
rz(-2.1814929) q[3];
sx q[3];
rz(-2.0898537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272098) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(0.89624727) q[0];
rz(-2.8889636) q[1];
sx q[1];
rz(-2.2815506) q[1];
sx q[1];
rz(-1.0505229) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79118585) q[0];
sx q[0];
rz(-1.5205654) q[0];
sx q[0];
rz(-2.6149261) q[0];
x q[1];
rz(-2.975198) q[2];
sx q[2];
rz(-0.67182589) q[2];
sx q[2];
rz(0.8478129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3747762) q[1];
sx q[1];
rz(-1.720806) q[1];
sx q[1];
rz(1.1039724) q[1];
x q[2];
rz(1.6231617) q[3];
sx q[3];
rz(-0.81365055) q[3];
sx q[3];
rz(0.46458515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.62319055) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(1.5790348) q[2];
rz(-1.7744428) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(2.3794543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52687454) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(2.7516464) q[0];
rz(0.82603106) q[1];
sx q[1];
rz(-2.4470058) q[1];
sx q[1];
rz(2.0894076) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74653502) q[0];
sx q[0];
rz(-1.87687) q[0];
sx q[0];
rz(2.2827882) q[0];
rz(-1.3953392) q[2];
sx q[2];
rz(-1.1793841) q[2];
sx q[2];
rz(3.1188426) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6893951) q[1];
sx q[1];
rz(-0.493825) q[1];
sx q[1];
rz(0.92269519) q[1];
rz(-2.1972606) q[3];
sx q[3];
rz(-1.9667224) q[3];
sx q[3];
rz(-2.348071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0803926) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(2.6123987) q[2];
rz(-0.75096327) q[3];
sx q[3];
rz(-0.96143985) q[3];
sx q[3];
rz(0.19415893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26218364) q[0];
sx q[0];
rz(-0.59469596) q[0];
sx q[0];
rz(-1.9770812) q[0];
rz(-1.2871845) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(-0.50416344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2874448) q[0];
sx q[0];
rz(-1.3519671) q[0];
sx q[0];
rz(1.3059421) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5515366) q[2];
sx q[2];
rz(-1.9474721) q[2];
sx q[2];
rz(-2.6835359) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.049147) q[1];
sx q[1];
rz(-1.7610252) q[1];
sx q[1];
rz(-1.8407497) q[1];
rz(0.88621288) q[3];
sx q[3];
rz(-1.4014114) q[3];
sx q[3];
rz(0.90506314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21504271) q[2];
sx q[2];
rz(-1.5554917) q[2];
sx q[2];
rz(3.1070993) q[2];
rz(1.5283594) q[3];
sx q[3];
rz(-0.72051636) q[3];
sx q[3];
rz(0.31149402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16019776) q[0];
sx q[0];
rz(-1.5481411) q[0];
sx q[0];
rz(2.7465469) q[0];
rz(-1.002671) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(-2.6832081) q[2];
sx q[2];
rz(-1.3946563) q[2];
sx q[2];
rz(-1.3763225) q[2];
rz(-0.18319753) q[3];
sx q[3];
rz(-2.8092794) q[3];
sx q[3];
rz(2.487779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
