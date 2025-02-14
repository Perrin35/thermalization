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
rz(0.74547493) q[0];
sx q[0];
rz(-0.59362721) q[0];
sx q[0];
rz(-2.797085) q[0];
rz(-1.966882) q[1];
sx q[1];
rz(-0.36829683) q[1];
sx q[1];
rz(-0.60454291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0812735) q[0];
sx q[0];
rz(-0.7501157) q[0];
sx q[0];
rz(1.4846205) q[0];
x q[1];
rz(-0.069434631) q[2];
sx q[2];
rz(-1.3926818) q[2];
sx q[2];
rz(1.3645862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82935059) q[1];
sx q[1];
rz(-1.568746) q[1];
sx q[1];
rz(-2.7228628) q[1];
x q[2];
rz(-0.80896583) q[3];
sx q[3];
rz(-0.74796092) q[3];
sx q[3];
rz(-0.35424074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1623666) q[2];
sx q[2];
rz(-1.484885) q[2];
sx q[2];
rz(-1.4244728) q[2];
rz(-0.10719565) q[3];
sx q[3];
rz(-0.19459477) q[3];
sx q[3];
rz(2.181982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2694038) q[0];
sx q[0];
rz(-0.93282455) q[0];
sx q[0];
rz(1.8316356) q[0];
rz(0.45920363) q[1];
sx q[1];
rz(-1.867086) q[1];
sx q[1];
rz(-0.31712636) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8882282) q[0];
sx q[0];
rz(-1.7087052) q[0];
sx q[0];
rz(-1.6186369) q[0];
x q[1];
rz(-0.049811157) q[2];
sx q[2];
rz(-1.5731817) q[2];
sx q[2];
rz(-2.002169) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9307946) q[1];
sx q[1];
rz(-0.55131492) q[1];
sx q[1];
rz(3.0153794) q[1];
rz(-pi) q[2];
rz(1.7789654) q[3];
sx q[3];
rz(-2.31593) q[3];
sx q[3];
rz(2.7221808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4294879) q[2];
sx q[2];
rz(-1.8450582) q[2];
sx q[2];
rz(-0.8075766) q[2];
rz(2.5782222) q[3];
sx q[3];
rz(-0.49401504) q[3];
sx q[3];
rz(2.0886683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5092369) q[0];
sx q[0];
rz(-2.95166) q[0];
sx q[0];
rz(-2.307039) q[0];
rz(-0.50476152) q[1];
sx q[1];
rz(-0.36634645) q[1];
sx q[1];
rz(-1.4588446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6989649) q[0];
sx q[0];
rz(-2.1219398) q[0];
sx q[0];
rz(-1.8056662) q[0];
x q[1];
rz(2.4065532) q[2];
sx q[2];
rz(-2.1872518) q[2];
sx q[2];
rz(2.7360345) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31400934) q[1];
sx q[1];
rz(-1.4928668) q[1];
sx q[1];
rz(-2.4185836) q[1];
x q[2];
rz(0.57679983) q[3];
sx q[3];
rz(-2.24772) q[3];
sx q[3];
rz(0.13595447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10996455) q[2];
sx q[2];
rz(-1.4723023) q[2];
sx q[2];
rz(2.7785981) q[2];
rz(1.1686769) q[3];
sx q[3];
rz(-0.76255885) q[3];
sx q[3];
rz(0.76538435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4204243) q[0];
sx q[0];
rz(-2.6154501) q[0];
sx q[0];
rz(2.5508733) q[0];
rz(-2.4144454) q[1];
sx q[1];
rz(-0.37494451) q[1];
sx q[1];
rz(0.66853833) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.062807) q[0];
sx q[0];
rz(-1.6123839) q[0];
sx q[0];
rz(1.6522626) q[0];
rz(2.3756251) q[2];
sx q[2];
rz(-1.9231638) q[2];
sx q[2];
rz(-2.3841928) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6741028) q[1];
sx q[1];
rz(-2.2558442) q[1];
sx q[1];
rz(-0.62226198) q[1];
x q[2];
rz(2.4005684) q[3];
sx q[3];
rz(-1.9771358) q[3];
sx q[3];
rz(2.5056553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.10355243) q[2];
sx q[2];
rz(-1.870564) q[2];
sx q[2];
rz(-1.6345778) q[2];
rz(-0.42190894) q[3];
sx q[3];
rz(-2.3814337) q[3];
sx q[3];
rz(-2.4222477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5108532) q[0];
sx q[0];
rz(-0.35950867) q[0];
sx q[0];
rz(1.0605633) q[0];
rz(0.063449055) q[1];
sx q[1];
rz(-2.5109406) q[1];
sx q[1];
rz(0.55387703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6543401) q[0];
sx q[0];
rz(-2.3759288) q[0];
sx q[0];
rz(-2.2310217) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2138344) q[2];
sx q[2];
rz(-1.0500488) q[2];
sx q[2];
rz(2.3152604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7859232) q[1];
sx q[1];
rz(-2.8357949) q[1];
sx q[1];
rz(3.0088916) q[1];
rz(-1.8964256) q[3];
sx q[3];
rz(-1.4417329) q[3];
sx q[3];
rz(-1.3323402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.267103) q[2];
sx q[2];
rz(-2.5449982) q[2];
sx q[2];
rz(-2.5388517) q[2];
rz(0.069325773) q[3];
sx q[3];
rz(-1.5963138) q[3];
sx q[3];
rz(-2.9749405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9423187) q[0];
sx q[0];
rz(-0.60180226) q[0];
sx q[0];
rz(0.45403516) q[0];
rz(-1.8858689) q[1];
sx q[1];
rz(-2.1818826) q[1];
sx q[1];
rz(1.7032636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3011264) q[0];
sx q[0];
rz(-2.2181803) q[0];
sx q[0];
rz(-2.6324138) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27737646) q[2];
sx q[2];
rz(-1.2834594) q[2];
sx q[2];
rz(-2.4871021) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2809873) q[1];
sx q[1];
rz(-1.0088723) q[1];
sx q[1];
rz(2.5562296) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71261974) q[3];
sx q[3];
rz(-1.1776867) q[3];
sx q[3];
rz(-2.676323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.85256514) q[2];
sx q[2];
rz(-1.270741) q[2];
sx q[2];
rz(1.7860058) q[2];
rz(1.9224613) q[3];
sx q[3];
rz(-1.4798603) q[3];
sx q[3];
rz(1.5662955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.3179625) q[0];
sx q[0];
rz(-0.83331236) q[0];
sx q[0];
rz(1.3939567) q[0];
rz(-0.45669237) q[1];
sx q[1];
rz(-0.87867457) q[1];
sx q[1];
rz(-1.3880233) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0972892) q[0];
sx q[0];
rz(-2.7345022) q[0];
sx q[0];
rz(1.2319698) q[0];
rz(-pi) q[1];
rz(-2.6235103) q[2];
sx q[2];
rz(-2.9112795) q[2];
sx q[2];
rz(-1.3467333) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2051831) q[1];
sx q[1];
rz(-0.35297063) q[1];
sx q[1];
rz(-2.754753) q[1];
x q[2];
rz(-2.0783992) q[3];
sx q[3];
rz(-1.0143544) q[3];
sx q[3];
rz(-2.0991652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2726511) q[2];
sx q[2];
rz(-0.77241263) q[2];
sx q[2];
rz(-0.19115494) q[2];
rz(-1.8027421) q[3];
sx q[3];
rz(-0.50722417) q[3];
sx q[3];
rz(0.52201456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1538447) q[0];
sx q[0];
rz(-2.4697883) q[0];
sx q[0];
rz(-1.9195358) q[0];
rz(-2.1568495) q[1];
sx q[1];
rz(-0.9877111) q[1];
sx q[1];
rz(-0.15880671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5903871) q[0];
sx q[0];
rz(-2.231601) q[0];
sx q[0];
rz(-1.4676276) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0350948) q[2];
sx q[2];
rz(-2.7125008) q[2];
sx q[2];
rz(0.96218357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8196311) q[1];
sx q[1];
rz(-2.0109977) q[1];
sx q[1];
rz(-2.0318713) q[1];
rz(1.430159) q[3];
sx q[3];
rz(-1.780073) q[3];
sx q[3];
rz(-1.9491553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1211991) q[2];
sx q[2];
rz(-2.8195916) q[2];
sx q[2];
rz(-0.89312345) q[2];
rz(-2.761306) q[3];
sx q[3];
rz(-1.2023353) q[3];
sx q[3];
rz(1.4343542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073008386) q[0];
sx q[0];
rz(-2.6717654) q[0];
sx q[0];
rz(1.481886) q[0];
rz(2.1549554) q[1];
sx q[1];
rz(-1.4886798) q[1];
sx q[1];
rz(2.5843487) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6300723) q[0];
sx q[0];
rz(-1.6247388) q[0];
sx q[0];
rz(0.11210895) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2115914) q[2];
sx q[2];
rz(-1.85382) q[2];
sx q[2];
rz(0.67365269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6861048) q[1];
sx q[1];
rz(-2.0596658) q[1];
sx q[1];
rz(-2.1829851) q[1];
x q[2];
rz(-2.557665) q[3];
sx q[3];
rz(-3.1213396) q[3];
sx q[3];
rz(-0.7113061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3096932) q[2];
sx q[2];
rz(-0.57277402) q[2];
sx q[2];
rz(2.0834303) q[2];
rz(1.7582827) q[3];
sx q[3];
rz(-1.0388831) q[3];
sx q[3];
rz(-0.96341187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-2.6365373) q[0];
sx q[0];
rz(-1.9872682) q[0];
sx q[0];
rz(-2.7099047) q[0];
rz(0.70026669) q[1];
sx q[1];
rz(-1.2338748) q[1];
sx q[1];
rz(-0.69469992) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0330677) q[0];
sx q[0];
rz(-1.0673079) q[0];
sx q[0];
rz(2.7208072) q[0];
rz(0.1044214) q[2];
sx q[2];
rz(-1.1200783) q[2];
sx q[2];
rz(1.3018276) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.34040117) q[1];
sx q[1];
rz(-0.71497289) q[1];
sx q[1];
rz(0.84960501) q[1];
rz(-pi) q[2];
rz(1.6142817) q[3];
sx q[3];
rz(-1.7491893) q[3];
sx q[3];
rz(2.1977193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6557287) q[2];
sx q[2];
rz(-2.8317917) q[2];
sx q[2];
rz(-2.7244205) q[2];
rz(-0.74472767) q[3];
sx q[3];
rz(-1.797902) q[3];
sx q[3];
rz(-0.97851306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3485296) q[0];
sx q[0];
rz(-1.2611669) q[0];
sx q[0];
rz(1.4766759) q[0];
rz(1.2186125) q[1];
sx q[1];
rz(-1.1398133) q[1];
sx q[1];
rz(-2.555991) q[1];
rz(2.2300874) q[2];
sx q[2];
rz(-2.3471033) q[2];
sx q[2];
rz(-2.4056619) q[2];
rz(-2.0176061) q[3];
sx q[3];
rz(-0.67025718) q[3];
sx q[3];
rz(2.0197301) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
