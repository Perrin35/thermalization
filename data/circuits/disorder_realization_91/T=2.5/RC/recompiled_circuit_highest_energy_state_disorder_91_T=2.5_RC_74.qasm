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
rz(-2.7741127) q[0];
sx q[0];
rz(-1.8125266) q[0];
sx q[0];
rz(1.654806) q[0];
rz(-1.6410671) q[1];
sx q[1];
rz(2.5379116) q[1];
sx q[1];
rz(11.707468) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6971815) q[0];
sx q[0];
rz(-1.0134122) q[0];
sx q[0];
rz(-2.0472705) q[0];
rz(-pi) q[1];
rz(0.5406981) q[2];
sx q[2];
rz(-1.0519111) q[2];
sx q[2];
rz(-0.93604871) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.20330435) q[1];
sx q[1];
rz(-1.4611725) q[1];
sx q[1];
rz(-2.2548864) q[1];
x q[2];
rz(-0.96894354) q[3];
sx q[3];
rz(-2.7730983) q[3];
sx q[3];
rz(0.70337765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19771244) q[2];
sx q[2];
rz(-0.79215017) q[2];
sx q[2];
rz(1.0666749) q[2];
rz(-0.094430447) q[3];
sx q[3];
rz(-2.5240199) q[3];
sx q[3];
rz(-0.42816952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18148947) q[0];
sx q[0];
rz(-0.26679376) q[0];
sx q[0];
rz(1.2130523) q[0];
rz(-0.45252291) q[1];
sx q[1];
rz(-2.8892543) q[1];
sx q[1];
rz(-2.4620893) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54645854) q[0];
sx q[0];
rz(-1.6041557) q[0];
sx q[0];
rz(0.24947687) q[0];
rz(-pi) q[1];
x q[1];
rz(1.405473) q[2];
sx q[2];
rz(-1.2836873) q[2];
sx q[2];
rz(1.6197255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7391239) q[1];
sx q[1];
rz(-1.1796412) q[1];
sx q[1];
rz(0.7372589) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4538058) q[3];
sx q[3];
rz(-1.7264328) q[3];
sx q[3];
rz(-2.9693672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54000336) q[2];
sx q[2];
rz(-0.82401472) q[2];
sx q[2];
rz(0.69671112) q[2];
rz(1.2434897) q[3];
sx q[3];
rz(-1.9466629) q[3];
sx q[3];
rz(-2.5534326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8884647) q[0];
sx q[0];
rz(-2.1051814) q[0];
sx q[0];
rz(-2.3179407) q[0];
rz(-0.15692391) q[1];
sx q[1];
rz(-1.7053968) q[1];
sx q[1];
rz(0.40208152) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21923301) q[0];
sx q[0];
rz(-1.439221) q[0];
sx q[0];
rz(1.1716362) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7355315) q[2];
sx q[2];
rz(-0.28898063) q[2];
sx q[2];
rz(-1.8584205) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3797219) q[1];
sx q[1];
rz(-2.5609697) q[1];
sx q[1];
rz(-2.1942744) q[1];
x q[2];
rz(-1.9605432) q[3];
sx q[3];
rz(-1.7754835) q[3];
sx q[3];
rz(1.2677416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0154401) q[2];
sx q[2];
rz(-0.93706477) q[2];
sx q[2];
rz(-0.13119571) q[2];
rz(-1.0570071) q[3];
sx q[3];
rz(-1.9481235) q[3];
sx q[3];
rz(-1.2098275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4256725) q[0];
sx q[0];
rz(-1.1765867) q[0];
sx q[0];
rz(-1.7536989) q[0];
rz(-1.7790986) q[1];
sx q[1];
rz(-1.8905996) q[1];
sx q[1];
rz(-1.2066134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20509991) q[0];
sx q[0];
rz(-1.1171891) q[0];
sx q[0];
rz(0.68286796) q[0];
rz(-pi) q[1];
rz(2.488606) q[2];
sx q[2];
rz(-0.76344244) q[2];
sx q[2];
rz(0.12432822) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.083664069) q[1];
sx q[1];
rz(-1.6477723) q[1];
sx q[1];
rz(-0.30268354) q[1];
x q[2];
rz(-0.73829262) q[3];
sx q[3];
rz(-1.2311934) q[3];
sx q[3];
rz(0.58635724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6917307) q[2];
sx q[2];
rz(-1.8064156) q[2];
sx q[2];
rz(0.94592363) q[2];
rz(-2.6876884) q[3];
sx q[3];
rz(-2.4920521) q[3];
sx q[3];
rz(-0.18979931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396963) q[0];
sx q[0];
rz(-1.8296158) q[0];
sx q[0];
rz(2.358118) q[0];
rz(-2.2354194) q[1];
sx q[1];
rz(-2.2515191) q[1];
sx q[1];
rz(-2.5772212) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63097341) q[0];
sx q[0];
rz(-1.6897908) q[0];
sx q[0];
rz(1.7158767) q[0];
rz(-pi) q[1];
rz(1.4323813) q[2];
sx q[2];
rz(-0.40988806) q[2];
sx q[2];
rz(2.5382535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5301492) q[1];
sx q[1];
rz(-2.2112729) q[1];
sx q[1];
rz(-2.3534858) q[1];
rz(-0.71471148) q[3];
sx q[3];
rz(-1.6922975) q[3];
sx q[3];
rz(2.7845259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.924661) q[2];
sx q[2];
rz(-1.4035808) q[2];
sx q[2];
rz(-0.04436392) q[2];
rz(-1.7313322) q[3];
sx q[3];
rz(-0.20874617) q[3];
sx q[3];
rz(1.9467719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57589543) q[0];
sx q[0];
rz(-0.79224753) q[0];
sx q[0];
rz(-0.76882291) q[0];
rz(2.5104751) q[1];
sx q[1];
rz(-1.9335856) q[1];
sx q[1];
rz(0.15359503) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0781411) q[0];
sx q[0];
rz(-0.65839863) q[0];
sx q[0];
rz(-1.0195579) q[0];
rz(-pi) q[1];
rz(-1.2943489) q[2];
sx q[2];
rz(-0.9381367) q[2];
sx q[2];
rz(-2.0273493) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5271618) q[1];
sx q[1];
rz(-1.2896184) q[1];
sx q[1];
rz(-1.6826732) q[1];
rz(0.098922313) q[3];
sx q[3];
rz(-0.77380816) q[3];
sx q[3];
rz(-0.1734002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.067387335) q[2];
sx q[2];
rz(-1.5978483) q[2];
sx q[2];
rz(0.81573168) q[2];
rz(-3.0873114) q[3];
sx q[3];
rz(-1.7547601) q[3];
sx q[3];
rz(1.774971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5197007) q[0];
sx q[0];
rz(-2.0760355) q[0];
sx q[0];
rz(0.52072293) q[0];
rz(2.0479274) q[1];
sx q[1];
rz(-2.2126074) q[1];
sx q[1];
rz(-1.3826133) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0313697) q[0];
sx q[0];
rz(-2.5532604) q[0];
sx q[0];
rz(-1.3084685) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6925631) q[2];
sx q[2];
rz(-0.64856883) q[2];
sx q[2];
rz(-2.8798333) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0573314) q[1];
sx q[1];
rz(-1.1597301) q[1];
sx q[1];
rz(0.89849679) q[1];
x q[2];
rz(-2.6723271) q[3];
sx q[3];
rz(-1.3949782) q[3];
sx q[3];
rz(-1.3202536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9426721) q[2];
sx q[2];
rz(-0.3519381) q[2];
sx q[2];
rz(-1.675763) q[2];
rz(-2.8525225) q[3];
sx q[3];
rz(-1.3705658) q[3];
sx q[3];
rz(1.2808965) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0220303) q[0];
sx q[0];
rz(-2.1896095) q[0];
sx q[0];
rz(-2.6493678) q[0];
rz(0.64552632) q[1];
sx q[1];
rz(-1.5461092) q[1];
sx q[1];
rz(2.2135977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4700025) q[0];
sx q[0];
rz(-0.3151463) q[0];
sx q[0];
rz(2.0285839) q[0];
rz(-pi) q[1];
rz(-0.3587175) q[2];
sx q[2];
rz(-1.9374829) q[2];
sx q[2];
rz(-2.3334954) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4175083) q[1];
sx q[1];
rz(-2.2615396) q[1];
sx q[1];
rz(-2.7907985) q[1];
x q[2];
rz(-1.2008529) q[3];
sx q[3];
rz(-1.2725012) q[3];
sx q[3];
rz(-1.4972794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8202028) q[2];
sx q[2];
rz(-1.6150183) q[2];
sx q[2];
rz(0.75902933) q[2];
rz(-1.2718893) q[3];
sx q[3];
rz(-1.3262409) q[3];
sx q[3];
rz(-0.71182865) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3439381) q[0];
sx q[0];
rz(-2.5929218) q[0];
sx q[0];
rz(-0.77242533) q[0];
rz(1.3683569) q[1];
sx q[1];
rz(-1.1027579) q[1];
sx q[1];
rz(2.8372884) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6876831) q[0];
sx q[0];
rz(-1.7237614) q[0];
sx q[0];
rz(1.3970333) q[0];
rz(-pi) q[1];
rz(1.3450121) q[2];
sx q[2];
rz(-1.567237) q[2];
sx q[2];
rz(0.6706711) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9615677) q[1];
sx q[1];
rz(-1.7503382) q[1];
sx q[1];
rz(2.3984329) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7469962) q[3];
sx q[3];
rz(-1.6871042) q[3];
sx q[3];
rz(1.968984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.407939) q[2];
sx q[2];
rz(-0.38241461) q[2];
sx q[2];
rz(1.1747053) q[2];
rz(2.5041653) q[3];
sx q[3];
rz(-1.2624242) q[3];
sx q[3];
rz(1.2229961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680854) q[0];
sx q[0];
rz(-2.8097476) q[0];
sx q[0];
rz(-0.77274957) q[0];
rz(-2.6840774) q[1];
sx q[1];
rz(-1.3950149) q[1];
sx q[1];
rz(1.4332019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49559418) q[0];
sx q[0];
rz(-0.80627215) q[0];
sx q[0];
rz(0.92732556) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.068763305) q[2];
sx q[2];
rz(-1.2085074) q[2];
sx q[2];
rz(-0.96093528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5262774) q[1];
sx q[1];
rz(-1.8377536) q[1];
sx q[1];
rz(-1.1039614) q[1];
x q[2];
rz(-2.5664342) q[3];
sx q[3];
rz(-1.0527826) q[3];
sx q[3];
rz(-0.94747551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2210803) q[2];
sx q[2];
rz(-2.1996193) q[2];
sx q[2];
rz(-0.67678779) q[2];
rz(-1.6109198) q[3];
sx q[3];
rz(-0.30083209) q[3];
sx q[3];
rz(1.0845186) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9250225) q[0];
sx q[0];
rz(-1.8393479) q[0];
sx q[0];
rz(1.7787697) q[0];
rz(-0.88381797) q[1];
sx q[1];
rz(-0.73278905) q[1];
sx q[1];
rz(-0.97272452) q[1];
rz(0.45223805) q[2];
sx q[2];
rz(-1.0483208) q[2];
sx q[2];
rz(2.0365289) q[2];
rz(2.6780405) q[3];
sx q[3];
rz(-1.3561854) q[3];
sx q[3];
rz(-0.88903313) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
