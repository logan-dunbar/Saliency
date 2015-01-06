function [stressLevel, timeLeft, out] = runLogansWorkDay(stressLevel, timeLeft)

	originalTimeLeft = timeLeft;
	while (timeLeft > 0 && (originalTimeLeft-timeLeft) < 24*60)
		timeLeft = timeLeft - 1;
		stressLevel = stressLevel + (originalTimeLeft/timeLeft)^(originalTimeLeft-timeLeft);
        out(originalTimeLeft-timeLeft) = stressLevel;
	end

end