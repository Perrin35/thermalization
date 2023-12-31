OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30259351) q[0];
sx q[0];
rz(-1.6380881) q[0];
sx q[0];
rz(1.5124613) q[0];
rz(-pi) q[1];
rz(1.2201266) q[2];
sx q[2];
rz(-1.7979243) q[2];
sx q[2];
rz(0.35958689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4343623) q[1];
sx q[1];
rz(-1.1303567) q[1];
sx q[1];
rz(2.4643154) q[1];
rz(2.1360374) q[3];
sx q[3];
rz(-1.894265) q[3];
sx q[3];
rz(0.934787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1258939) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(-1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.9803479) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.6764486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11008308) q[0];
sx q[0];
rz(-2.2950036) q[0];
sx q[0];
rz(-1.5741041) q[0];
rz(0.25603489) q[2];
sx q[2];
rz(-1.0962152) q[2];
sx q[2];
rz(1.408996) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2966753) q[1];
sx q[1];
rz(-2.1293318) q[1];
sx q[1];
rz(-1.0122453) q[1];
x q[2];
rz(2.7705455) q[3];
sx q[3];
rz(-2.3608532) q[3];
sx q[3];
rz(2.3649529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(3.0318276) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(-1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(0.51613581) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-0.80054545) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0101937) q[0];
sx q[0];
rz(-2.7079765) q[0];
sx q[0];
rz(1.9171159) q[0];
x q[1];
rz(0.86513743) q[2];
sx q[2];
rz(-1.0312928) q[2];
sx q[2];
rz(1.3441966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86299455) q[1];
sx q[1];
rz(-1.1578214) q[1];
sx q[1];
rz(1.5763271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9242026) q[3];
sx q[3];
rz(-0.73261515) q[3];
sx q[3];
rz(0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(1.3595954) q[2];
rz(0.96757403) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(-2.676679) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(-2.0565313) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5497919) q[0];
sx q[0];
rz(-1.340197) q[0];
sx q[0];
rz(-2.6016298) q[0];
rz(-pi) q[1];
rz(-2.8934946) q[2];
sx q[2];
rz(-1.2144226) q[2];
sx q[2];
rz(-0.46596371) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9492053) q[1];
sx q[1];
rz(-2.9678223) q[1];
sx q[1];
rz(-2.5237571) q[1];
x q[2];
rz(-1.163108) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(2.98416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1365635) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(2.1599105) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0504793) q[0];
sx q[0];
rz(-2.1149153) q[0];
sx q[0];
rz(-2.1168461) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1443411) q[2];
sx q[2];
rz(-0.28595668) q[2];
sx q[2];
rz(-2.2821102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.423646) q[1];
sx q[1];
rz(-1.4232993) q[1];
sx q[1];
rz(2.6731078) q[1];
rz(-pi) q[2];
rz(-0.97335191) q[3];
sx q[3];
rz(-1.5095599) q[3];
sx q[3];
rz(2.6254568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(2.6521818) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(1.4846102) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(0.16528027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76089345) q[0];
sx q[0];
rz(-1.5415511) q[0];
sx q[0];
rz(0.0951035) q[0];
x q[1];
rz(-3.0566453) q[2];
sx q[2];
rz(-2.4521378) q[2];
sx q[2];
rz(1.2598318) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89989923) q[1];
sx q[1];
rz(-0.69679931) q[1];
sx q[1];
rz(-0.65710575) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0201488) q[3];
sx q[3];
rz(-1.7835788) q[3];
sx q[3];
rz(-1.3495812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2281987) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.4432663) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3570324) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(0.7473942) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-0.73648891) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38024494) q[0];
sx q[0];
rz(-1.5878829) q[0];
sx q[0];
rz(-1.6197617) q[0];
rz(2.1963504) q[2];
sx q[2];
rz(-0.40996273) q[2];
sx q[2];
rz(-1.9117102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10662096) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(-0.76645318) q[1];
rz(-1.1398846) q[3];
sx q[3];
rz(-1.1093372) q[3];
sx q[3];
rz(1.0678837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0050469) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(2.6100256) q[2];
rz(-0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(-1.846107) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078995973) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.2980365) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(0.92179006) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73496504) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(-0.69358967) q[0];
rz(-pi) q[1];
rz(1.8838896) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(1.8956172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7505155) q[1];
sx q[1];
rz(-0.32954307) q[1];
sx q[1];
rz(-2.4604172) q[1];
rz(0.86254085) q[3];
sx q[3];
rz(-1.7077912) q[3];
sx q[3];
rz(2.9555637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90199295) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.1716589) q[2];
rz(1.3575859) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(0.67970651) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4363791) q[0];
sx q[0];
rz(-2.0853015) q[0];
sx q[0];
rz(1.8041457) q[0];
x q[1];
rz(0.87551261) q[2];
sx q[2];
rz(-1.4317703) q[2];
sx q[2];
rz(-1.306844) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19787994) q[1];
sx q[1];
rz(-2.3730952) q[1];
sx q[1];
rz(-0.49853034) q[1];
rz(-2.4217442) q[3];
sx q[3];
rz(-0.82092972) q[3];
sx q[3];
rz(-1.6380701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1264964) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(2.2311907) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0891721) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(-1.9649327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0929716) q[0];
sx q[0];
rz(-1.7341359) q[0];
sx q[0];
rz(-3.0924348) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2048336) q[2];
sx q[2];
rz(-2.0271795) q[2];
sx q[2];
rz(2.2697743) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9709819) q[1];
sx q[1];
rz(-0.5628399) q[1];
sx q[1];
rz(1.7112205) q[1];
rz(-pi) q[2];
rz(0.51104607) q[3];
sx q[3];
rz(-1.647445) q[3];
sx q[3];
rz(-2.6443036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(-2.0142377) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-2.4049315) q[2];
sx q[2];
rz(-1.7792637) q[2];
sx q[2];
rz(-1.7532495) q[2];
rz(0.022396537) q[3];
sx q[3];
rz(-0.35084421) q[3];
sx q[3];
rz(2.4435333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
