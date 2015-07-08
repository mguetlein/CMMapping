package org.chesmapper.map.alg;

public class RuntimeOwner
{
	private long runtimeStart;
	private long runtime;

	protected void startRuntime()
	{
		runtimeStart = System.currentTimeMillis();
	}

	protected void stopRuntime()
	{
		if (runtimeStart == 0)
			throw new IllegalStateException();
		runtime = System.currentTimeMillis() - runtimeStart;
	}

	public long getRuntime()
	{
		return runtime;
	}
}
