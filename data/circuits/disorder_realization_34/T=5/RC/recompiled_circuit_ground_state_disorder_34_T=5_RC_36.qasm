OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4021969) q[0];
sx q[0];
rz(-0.88391179) q[0];
sx q[0];
rz(0.85564268) q[0];
rz(-0.17172509) q[1];
sx q[1];
rz(6.167616) q[1];
sx q[1];
rz(12.003916) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2320975) q[0];
sx q[0];
rz(-1.9350855) q[0];
sx q[0];
rz(-3.0845736) q[0];
x q[1];
rz(1.2954584) q[2];
sx q[2];
rz(-0.52342969) q[2];
sx q[2];
rz(1.1688423) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89002883) q[1];
sx q[1];
rz(-1.152413) q[1];
sx q[1];
rz(-2.0189925) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3381932) q[3];
sx q[3];
rz(-0.53871846) q[3];
sx q[3];
rz(2.5634457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2354551) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(-2.3423024) q[2];
rz(-0.47131395) q[3];
sx q[3];
rz(-0.91336942) q[3];
sx q[3];
rz(-0.93588626) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1021295) q[0];
sx q[0];
rz(-1.3251745) q[0];
sx q[0];
rz(-1.7934196) q[0];
rz(-1.3735636) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(-1.9893533) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9236167) q[0];
sx q[0];
rz(-2.2647175) q[0];
sx q[0];
rz(-2.7569735) q[0];
x q[1];
rz(-2.197108) q[2];
sx q[2];
rz(-2.1462893) q[2];
sx q[2];
rz(1.4552417) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2272646) q[1];
sx q[1];
rz(-0.89114571) q[1];
sx q[1];
rz(1.3548681) q[1];
x q[2];
rz(2.3186734) q[3];
sx q[3];
rz(-2.4429818) q[3];
sx q[3];
rz(0.60958344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3422602) q[2];
sx q[2];
rz(-0.41877425) q[2];
sx q[2];
rz(0.93117923) q[2];
rz(3.0350507) q[3];
sx q[3];
rz(-2.0276766) q[3];
sx q[3];
rz(-2.54134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057864144) q[0];
sx q[0];
rz(-1.3902384) q[0];
sx q[0];
rz(2.6237543) q[0];
rz(-0.88874108) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(-2.6944366) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5260552) q[0];
sx q[0];
rz(-1.8319905) q[0];
sx q[0];
rz(-1.636807) q[0];
rz(2.1250399) q[2];
sx q[2];
rz(-2.2187761) q[2];
sx q[2];
rz(-2.3584443) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7501156) q[1];
sx q[1];
rz(-0.98408723) q[1];
sx q[1];
rz(2.9787885) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8497883) q[3];
sx q[3];
rz(-1.0615665) q[3];
sx q[3];
rz(-1.694547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37457028) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(2.6089597) q[2];
rz(2.8940708) q[3];
sx q[3];
rz(-0.73740021) q[3];
sx q[3];
rz(2.1029162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5175051) q[0];
sx q[0];
rz(-3.0563323) q[0];
sx q[0];
rz(3.0349773) q[0];
rz(2.7952349) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(1.1603629) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7713523) q[0];
sx q[0];
rz(-1.816907) q[0];
sx q[0];
rz(-1.3918124) q[0];
rz(-2.423942) q[2];
sx q[2];
rz(-2.1776878) q[2];
sx q[2];
rz(-0.90536149) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6667156) q[1];
sx q[1];
rz(-1.2704388) q[1];
sx q[1];
rz(-1.9431861) q[1];
rz(0.45063536) q[3];
sx q[3];
rz(-1.0900146) q[3];
sx q[3];
rz(-1.4522875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13566636) q[2];
sx q[2];
rz(-0.29834193) q[2];
sx q[2];
rz(3.0653817) q[2];
rz(2.5557319) q[3];
sx q[3];
rz(-1.9867089) q[3];
sx q[3];
rz(1.3425672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0214486) q[0];
sx q[0];
rz(-1.8360538) q[0];
sx q[0];
rz(-0.35010499) q[0];
rz(2.1856951) q[1];
sx q[1];
rz(-1.2799542) q[1];
sx q[1];
rz(-1.9116481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53287017) q[0];
sx q[0];
rz(-1.2048831) q[0];
sx q[0];
rz(-0.25718148) q[0];
rz(0.036418865) q[2];
sx q[2];
rz(-1.4055739) q[2];
sx q[2];
rz(-0.18319229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0240682) q[1];
sx q[1];
rz(-1.9000016) q[1];
sx q[1];
rz(-0.19964053) q[1];
x q[2];
rz(1.2754945) q[3];
sx q[3];
rz(-2.5341138) q[3];
sx q[3];
rz(-2.5584084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.89796394) q[2];
sx q[2];
rz(-0.25883365) q[2];
sx q[2];
rz(-1.6492856) q[2];
rz(-2.7367075) q[3];
sx q[3];
rz(-0.76571524) q[3];
sx q[3];
rz(-2.305472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355857) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(-0.58240044) q[0];
rz(2.4549585) q[1];
sx q[1];
rz(-1.4056987) q[1];
sx q[1];
rz(-1.883421) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67418881) q[0];
sx q[0];
rz(-1.4454525) q[0];
sx q[0];
rz(2.0095429) q[0];
x q[1];
rz(2.7863726) q[2];
sx q[2];
rz(-2.0783011) q[2];
sx q[2];
rz(-1.2025646) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2542779) q[1];
sx q[1];
rz(-1.5149024) q[1];
sx q[1];
rz(0.069737597) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9223676) q[3];
sx q[3];
rz(-2.1532183) q[3];
sx q[3];
rz(-1.6265914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76576343) q[2];
sx q[2];
rz(-1.148369) q[2];
sx q[2];
rz(-1.0180391) q[2];
rz(-2.8954519) q[3];
sx q[3];
rz(-1.3956416) q[3];
sx q[3];
rz(1.5181946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.70596424) q[0];
sx q[0];
rz(-0.80393296) q[0];
sx q[0];
rz(2.5614118) q[0];
rz(-2.9980581) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(2.9339583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.62791) q[0];
sx q[0];
rz(-1.915689) q[0];
sx q[0];
rz(0.78457997) q[0];
x q[1];
rz(-2.7933664) q[2];
sx q[2];
rz(-2.1219606) q[2];
sx q[2];
rz(-1.4081692) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.27397284) q[1];
sx q[1];
rz(-1.8302396) q[1];
sx q[1];
rz(-2.2745489) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36334857) q[3];
sx q[3];
rz(-1.4814113) q[3];
sx q[3];
rz(2.9575728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84270728) q[2];
sx q[2];
rz(-1.8536114) q[2];
sx q[2];
rz(2.5353954) q[2];
rz(2.4541564) q[3];
sx q[3];
rz(-2.0076553) q[3];
sx q[3];
rz(-2.4351951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48924482) q[0];
sx q[0];
rz(-0.027996538) q[0];
sx q[0];
rz(1.0580753) q[0];
rz(-0.10969133) q[1];
sx q[1];
rz(-1.1166162) q[1];
sx q[1];
rz(-1.6995957) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083695166) q[0];
sx q[0];
rz(-2.6287615) q[0];
sx q[0];
rz(-0.67016853) q[0];
x q[1];
rz(-0.47693129) q[2];
sx q[2];
rz(-2.5897346) q[2];
sx q[2];
rz(0.98699206) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86275872) q[1];
sx q[1];
rz(-2.2087084) q[1];
sx q[1];
rz(2.5274123) q[1];
rz(1.2992925) q[3];
sx q[3];
rz(-1.3437004) q[3];
sx q[3];
rz(0.61448594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.480392) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(1.2083758) q[2];
rz(2.4783573) q[3];
sx q[3];
rz(-1.4570718) q[3];
sx q[3];
rz(1.2164345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1709568) q[0];
sx q[0];
rz(-2.1747776) q[0];
sx q[0];
rz(3.048625) q[0];
rz(1.8611106) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(0.13519898) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4358035) q[0];
sx q[0];
rz(-3.0717141) q[0];
sx q[0];
rz(0.70894928) q[0];
rz(-2.8835758) q[2];
sx q[2];
rz(-1.7191559) q[2];
sx q[2];
rz(2.6515863) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0307903) q[1];
sx q[1];
rz(-0.25486481) q[1];
sx q[1];
rz(0.084666208) q[1];
rz(-1.7235123) q[3];
sx q[3];
rz(-1.6441321) q[3];
sx q[3];
rz(-0.45460816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70904237) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(-0.77152983) q[2];
rz(0.49154526) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(2.5468723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5563357) q[0];
sx q[0];
rz(-2.519643) q[0];
sx q[0];
rz(2.0167895) q[0];
rz(1.2987761) q[1];
sx q[1];
rz(-2.5262084) q[1];
sx q[1];
rz(0.76464701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52228868) q[0];
sx q[0];
rz(-1.4656656) q[0];
sx q[0];
rz(-1.6572857) q[0];
x q[1];
rz(-0.89057335) q[2];
sx q[2];
rz(-1.169551) q[2];
sx q[2];
rz(0.97637343) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54658871) q[1];
sx q[1];
rz(-1.6975132) q[1];
sx q[1];
rz(1.3664043) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9493136) q[3];
sx q[3];
rz(-1.0487742) q[3];
sx q[3];
rz(1.4498364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3760959) q[2];
sx q[2];
rz(-1.2906047) q[2];
sx q[2];
rz(-0.20720227) q[2];
rz(0.97992212) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(-2.9492212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021066396) q[0];
sx q[0];
rz(-1.4561894) q[0];
sx q[0];
rz(-0.86984632) q[0];
rz(-2.6407241) q[1];
sx q[1];
rz(-2.905838) q[1];
sx q[1];
rz(0.9314608) q[1];
rz(2.4293368) q[2];
sx q[2];
rz(-2.9604572) q[2];
sx q[2];
rz(-1.0618718) q[2];
rz(2.9034782) q[3];
sx q[3];
rz(-1.7287711) q[3];
sx q[3];
rz(0.14915376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
